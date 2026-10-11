# The dataset's version-3 texts: native observations before the reader takes their header

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0.
`compare_brep_io.py --capture-v3 --v3-dataset DIR` found the public
dataset's `.brep` texts whose header is OCCT's version-3 line (`CASCADE
Topology V3, (c) Open Cascade`, matched as `TopTools_ShapeSet::Read` matches
it) and read each with the two probes `compare_brep_io.py` already pins,
unchanged (`occt_brep_io_oracle.cpp` for solids, `occt_free_shape_oracle.cpp`
for free shapes; `capture.json` records both sources' digests). The dataset
is local-only (U1): this directory holds the files' relative paths and
digests and the native observations, never a file's text.

`capture.json` records the Rust revision and the kernel's reading of every
file before any code took the line (`brep_io_probe read`): each refused at
the header, `expected Matra-Datavision` (`rust_v3_header_read` false).
`native-solids.txt` lists, per file, every solid reached from the root
through compounds with its BRepCheck verdict, distinct subshape counts and
volume, area and centre; `native-free.txt` every shell, face, wire, edge and
vertex outside a solid with its type, verdict, counts and area or length
with its centre.

## Observations

* 26 texts (25 in `brep/`, `others/bug30175_box.brep`); OCCT reads every
  one.
* 11 solids, all BRepCheck-valid: `bug28620` (1), `bug32942` to `bug32945`
  (2 each) and `bug32960` (2).
* 28 free shapes: 20 faces, 6 wires and 2 edges. Valid but for 9 faces:
  `bug27170_f`'s, `bug32617`'s six (faces without a surface on
  triangulations with normals), `bug32737_free_verts`'s and
  `bug30175_box`'s (a face without a surface or edges on a triangulation;
  BRepGProp measures its triangles).
