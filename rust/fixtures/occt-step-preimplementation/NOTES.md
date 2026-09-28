# STEP files read natively (STEP-a): native observations before the importer

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks and the data exchange toolkits
(`build_pinned_occt.py --toolkit TKDESTEP`). `capture.json` records the Rust
revision and that the kernel's STEP importer did not exist
(`rust/kernel/src/step/import.rs` absent; the Part 21 reader
`step/part21.rs` and its tests did). `inputs.txt` holds the 22 fixture files
of `fixtures/step` (`generate_step_fixtures.py`), each framed as `case NAME
LENGTH` and its bytes; `oracle.cpp` is `occt_step_oracle.cpp` as captured:
each file read by `STEPControl_Reader::ReadStream` at the reader's defaults,
every root transferred, and `native.txt` its bodies: class, OCCT's counts,
BRepCheck's verdict, volume, area, centre and largest tolerance.

## Observations

* Every file reads, every body is valid, and every body is the independent
  reference's (`step-expected.tsv`): the same class and counts (a cone's apex
  and a sphere's poles each gain the degenerated edge the reader adds, the
  torus has one vertex and two seam edges), volumes, areas and centres
  within `7.6e-15` relative (the sphere; planar and cylindrical bodies agree
  to `2.5e-16`).
* Units: the box in metres and in inches reads in millimetres, the
  frustum's semi-angle in degrees as the one in radians; the hand-written
  file (comments, complex instances out of order, forward references,
  unnormalised directions, `$` reference directions, a face bound of
  orientation `.F.` on a surface of the opposite normal) reads as the
  corner tetrahedron.
* Every tolerance after the reader's shape healing is `1e-7` mm, the files'
  uncertainty: nothing was enlarged.

`compare_step.py` requires every later run to reproduce these rows (on
another platform, its reviewed record) and compares the kernel's bodies with
both the reference and these observations.
