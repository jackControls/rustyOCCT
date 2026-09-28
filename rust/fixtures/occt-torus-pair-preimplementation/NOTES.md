# Two tori off a common axis (S7b.3b.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that the kernel did not intersect two tori
(`rust/kernel/src/intersection/torus_curves.rs` without `torus_pair`).
`inputs.txt` holds the 14 `tt_` cases of `torus-curve-cases.txt` with the
kernel's stored frame axes; `oracle.cpp` is
`occt_procedural_intersection_oracle.cpp` as captured, and `native.txt` its
rows: `GeomInt_IntSS` with tolerance `1e-7`, each line's closedness and 17
samples, and each isolated point.

## Observations

* Against the independent reference (`torus_curve_reference.py` with the
  other torus's implicit quartic, whose restriction to a meridian circle is
  still of degree two in `cos t`, `sin t`; tangencies from a Groebner basis
  of the two spines' critical pairs), 10 cases agree: every native sample
  within `1e-6` of both tori, every component covered, the empty case empty.
* Reviewed (`occt-torus-curve-divergences.json`): three isolated tangencies
  give no native point (outer equators, a saddle touched from inside, a tube
  resting on another at right angles), and the loop of two tori `2^-20`
  closer than touching is missed.

`compare_torus_curves.py` requires every later run to reproduce these rows
(on another platform, its reviewed record) and compares the kernel's results
with both.
