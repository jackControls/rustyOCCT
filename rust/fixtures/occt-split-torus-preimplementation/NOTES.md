# Tori split by a plane (S8d): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for the split of tori existed
(`rust/kernel/src/solid/split/torus.rs` absent). `inputs.txt` holds the 12
cases of `split-torus-cases.txt`, each a whole torus made by
`BRepPrimAPI_MakeTorus` on a `gp_Ax2` with the kernel's stored frame axes,
with the splitting plane; `oracle.cpp` is `occt_split_oracle.cpp` as
captured, and `native.txt` its solids: side, volume, area, centre, face,
edge and vertex counts, validity.

## Observations

* Against the independent reference (`split_reference.torus_rows`: each
  slice an annulus cut by the plane's line), the cuts normal to the axis
  or containing it agree to 1e-15 relative. On the planes cutting the tube
  in spiric curves (parallel to the axis through the tube, or oblique),
  OCCT's section edges and faces are B-spline approximations and
  BRepGProp's volumes, areas and centres differ from the reference by up to
  3e-7 relative: reviewed native differences once the kernel splits them.
* Every piece is valid. A plane touching the tube's top or missing it
  leaves one solid.
* OCCT keeps the torus's seams (at `u = 0` and `v = 0`), so its counts
  differ from the kernel's seamless cells wherever a piece keeps one.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record).
