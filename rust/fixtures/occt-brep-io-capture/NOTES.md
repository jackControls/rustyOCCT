# .brep interop (T2): native observations of the upstream corpus

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0 built
headless with exception checks (`build_pinned_occt.py`, toolkits through
TKTopAlgo). `oracle.cpp` is `rust/tools/occt_brep_io_oracle.cpp` as captured;
`native.txt` holds its output for every `data/occ/*.brep` file, unmodified
from upstream (`capture.json` pins each file's SHA-256).

## Order of evidence: this capture came after the implementation

Every earlier milestone captured native observations before any Rust code
existed. T2 did not. The kernel's reader, converter and writer
(`rust/kernel/src/occt_brep*`) were written first and were developed
against the corpus files themselves, then this capture was taken on
2026-09-26 in the same working tree. The break is recorded here and not
repaired, because it cannot be repaired after the fact. What it does and
does not affect:

* The captured observations depend on no Rust output. They are OCCT
  reading unmodified upstream files: which solids the root reaches, their
  BRepCheck verdicts, distinct subshape counts and mass properties. Any
  later run must reproduce them (counts and verdicts exactly, properties
  within `PROPERTY_BOUND`).
* The expectations the kernel is tested against come from the independent
  reader `brep_io_reference.py`. It was written from the format
  specification (`dox/specification/brep_format.md`) and the pinned
  `BRepTools_ShapeSet` / `TopTools_LocationSet` sources, without Rust. It
  was certified against native OCCT only after both existed. It was not
  derived from the native observations, but it was written by the same
  author who had already debugged the Rust reader on these files, so it is
  not blind to what that work revealed (see below).
* The round trips are observed fresh on every run and never captured:
  every identity prism, and every corpus solid the kernel imports, as the
  kernel writes them. Their expectations come from the kernel's exact mass
  properties (prisms) and from this capture (corpus solids), not from the
  written text.

## What the corpus revealed about the format

Facts the Rust reader first got wrong on these files, all confirmed in the
pinned sources:

* A matrix location record (`1`) always takes the next index, even when it
  is the identity. An empty composite chain (`2 0`) is not numbered.
  Composition is `L = L1^p * L` in record order.
* Shape records are numbered backwards: the last record written is 1.
* A closed surface's continuity is written straight after the second
  pcurve number with no space (`3 12 13CN ...` becomes `13CN`), and
  `operator>>` reads it back as two tokens.
* B-spline curve and surface records carry rational and periodic flags
  before their degrees.
* UV end points on a pcurve representation appear only in format version 2
  exactly.
* A trimmed line or circle is its basis over a range; planes may omit
  pcurves (`BRep_Tool::CurveOnPlane` derives them).
* Cylinders can be indirect (`X ^ Y = -N`). Wires are unordered, and a
  face's outer wire is not necessarily first.
* Closed circles are written with range `6.28318530717959`, not `2 pi`, so
  pcurve ranges must be snapped before comparing them with `2 pi`.

## Coverage of this corpus

37 files, 77 solids reached from their roots. 29 solids have
structure the kernel represents (plane and cylinder faces, line and circle
edges, rigid locations). The kernel imports exactly those 29, and OCCT gives
their originals the counts the kernel synthesizes. Every other solid, and
every table record the kernel cannot represent, is reported by name and
count (`brep-io-expected.tsv`).
