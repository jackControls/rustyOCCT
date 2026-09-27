# The corpus's free shapes (S6): native observations before any kernel code imports them

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0.
`capture.json` records the Rust revision and that no kernel code imported a
shell, face, wire, edge or vertex outside a solid: the reader reported them
as `FreeShell`, `FreeFace`, `FreeWire`, `FreeEdge` and `FreeVertex`,
unsupported. `oracle.cpp` is `occt_free_shape_oracle.cpp` as captured (a
separate probe: the T2 capture pins `occt_brep_io_oracle.cpp`'s text), and
`native.txt` lists, per `data/occ` file, every free shape reached from the
root through compounds with its type, BRepCheck verdict, distinct subshape
counts, and area or length with its centre.

## Observations

* 6,223 free shapes in 12 files: 59 faces (`MAT`, `Room`, `face`, `face1`,
  `face2`, `hammer`, `terrain`, `wing`), 2 shells (`fuse`, `shell1`), 2,829
  edges (`edge`, and 2,828 in `asahi`) and 3,333 vertices (`asahi`). Every one
  is BRepCheck-valid.
* The independent reader (`brep_io_reference.free_shapes`) enumerates the
  same shapes with the same types and counts; `compare_brep_io.py
  --capture-free` refused to write the capture otherwise. By its rules
  `face2` (an ellipse pcurve on a trimmed surface) and `shell1` are not
  representable; every other free shape is.
