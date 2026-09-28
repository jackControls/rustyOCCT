# Cones and zones cut in conics (S8d.2): native observations after the kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. Unlike the other split captures this
one was taken after `rust/kernel/src/solid/split/conic.rs` existed, and
`capture.json` says so (`rust_split_exists: true`); `compare_split.py`
accepts it only under its post-implementation key (`s8d2`). The evidence
taken before any code for these sections is S8c's capture
(`occt-split-primitive-preimplementation`), which already held three of
them (`apex_oblique`, `frustum_parallel`, `zone_oblique`). These 12 cases of
`split-conic-cases.txt` widen the configurations (a closed ellipse, a
parabola, hyperbolas across one rim and both, a tongue at the top, the
rulings through a frustum's virtual apex, an ellipse touching both rims,
tilted cones and zones, a zone's side cap, a cap with its pole, circle arcs
across both rims); they were chosen from the configurations, not from the
kernel's results. `inputs.txt` holds them as made by `BRepPrimAPI_MakeCone`
and `MakeSphere` on the kernel's stored frame axes, `oracle.cpp` is
`occt_split_oracle.cpp` as captured and `native.txt` its solids.

## Observations

* Against the independent reference (`split_reference.revolved_rows`),
  every case agrees within 2e-8 but `cap_oblique`, whose spherical faces
  BRepGProp measures to 6.4e-8 relative (as `sphere_oblique`'s).
* Every piece is valid. OCCT keeps the seam at `u = 0` (the frame's x
  axis), splitting any face or edge of a piece it crosses; the kernel's
  seamless counts differ wherever it does (reviewed in
  `occt-split-divergences.json`).

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record).
