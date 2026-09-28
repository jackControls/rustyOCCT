# Cones and spheres split by a plane (S8c): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for the split of revolved solids existed
(`rust/kernel/src/solid/split/revolved.rs` absent). `inputs.txt` holds the 23
cases of `split-primitive-cases.txt`, each a `BRepPrimAPI_MakeCone` or
`BRepPrimAPI_MakeSphere` on a `gp_Ax2` with the kernel's stored frame axes
(`stored_axes`), with the splitting plane; `oracle.cpp` is
`occt_split_oracle.cpp` as captured, and `native.txt` its solids: side,
volume, area, centre, face, edge and vertex counts, validity.

## Observations

* Against the independent reference (`split_reference.revolved_rows`: each
  slice a disc cut by the plane's line, integrated along the axis), the
  cones, frusta and zones cut by planes normal to or containing the axis
  agree to 4e-15 relative. On cuts across the surface at an angle,
  BRepGProp's error grows: 4.2e-9 (`apex_oblique`), 1.3e-8
  (`zone_oblique`), 6.6e-8 (`sphere_oblique`), and 4.7e-4 on
  `sphere_tilted`, whose two pieces sum to 14.138748 where a sphere of
  radius 1.5 has 14.137167 (`4/3 pi 1.5^3`). Those are reviewed native
  differences, with the closed forms of a spherical cap as the independent
  evidence.
* Every piece is valid. A plane through an apex, lying in a cap, touching a
  sphere or missing leaves one solid.
* OCCT keeps each revolved face's seam at the frame's x axis (angle 0) and
  a degenerated edge at an apex or pole, so its counts differ from the
  kernel's seamless cells wherever a piece keeps a seam: reviewed per case
  once the kernel splits it.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record).
