# Tori cut in spiric sections (S8d.3): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for the spiric split existed
(`rust/kernel/src/solid/split/spiric.rs` absent). `inputs.txt` holds the 13
cases of `split-spiric-cases.txt`, each a whole torus made by
`BRepPrimAPI_MakeTorus` on a `gp_Ax2` with the kernel's stored frame axes,
with the splitting plane: two loops winding about the axis (`torus_gentle`,
`torus_band_offset`, `torus_top_band`, `torus_tilted_band`), one
contractible loop (`torus_top_cap`, `torus_cap_outer`, `torus_peanut`,
`torus_hole_slice`, `torus_tilted_spiric`) and two loops winding about the
tube (`torus_two_ovals`, `torus_steep`, `torus_skew_ovals`,
`torus_small_tube`), each checked by sampling the section's sign on the
torus before it was chosen. `oracle.cpp` is `occt_split_oracle.cpp` as
captured, and `native.txt` its solids.

## Observations

* Every piece is valid, and the sides agree with the reference
  (`split_reference.torus_rows`).
* OCCT's section edges and faces are B-spline approximations: BRepGProp's
  volumes and areas differ from the reference by 1.6e-7 to 9.2e-6 relative,
  and by 7.6e-4 on `torus_small_tube` (major 4, minor 0.5) and 6.9e-4 on
  `torus_tilted_spiric`: reviewed native differences once the kernel splits
  them.
* OCCT keeps the torus's seams (at `u = 0` and `v = 0`), so its counts
  differ from the kernel's seamless cells wherever a piece keeps one.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record).
