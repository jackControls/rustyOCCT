# Procedural intersections with cones (S7b.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for cone pairs existed
(`rust/kernel/src/intersection/procedural.rs` without `Surface::Cone`).
`inputs.txt` holds the 15 cone cases of `procedural-intersection-cases.txt`
(`ck_` a cylinder and a cone whose axes are not coaxial, `ks_` a sphere and a
cone with the centre off the axis) with the kernel's stored frame axes;
`oracle.cpp` is `occt_procedural_intersection_oracle.cpp` with cones (the
S7b.1 capture's probe, extended; its cylinder and sphere output is
unchanged), and `native.txt` its `GeomInt_IntSS` rows.

## Observations

* Against the independent reference, 14 cases agree: every native sample
  within `1e-6` of the case's size from the exact curve by closest-point
  distance (near a loop's end the branches are vertical in the ruled
  parameterisation, so a point's own angle is not its nearest parameter), and
  every component (each loop, each ring) carrying native samples.
* Reviewed (`occt-procedural-intersection-divergences.json`): a sphere
  `2^-20` beyond tangency to both nappes gives two tiny ill-conditioned
  loops, and one of OCCT's samples, within `2e-9` of both surfaces, is
  `1.55e-6` from the exact curve.
