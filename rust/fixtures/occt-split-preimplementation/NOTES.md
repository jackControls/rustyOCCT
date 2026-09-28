# Prisms split by a plane (S8a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks. `capture.json` records the Rust
revision and that no kernel code for the split existed
(`rust/kernel/src/solid/split.rs` absent). `inputs.txt` holds the 26 cases of
`split-cases.txt` as explicit constructions (`identity_reference.native_case`:
the profile face at the start offset, its prism vector) with the splitting
plane; `oracle.cpp` is `occt_split_oracle.cpp` as captured: the prism split
by `BRepAlgoAPI_Splitter` with a planar face 2,000 wide, and `native.txt` its
solids: side, volume, area, centre, face, edge and vertex counts, validity.

## Observations

* Against the independent reference (`split_reference.py`: each side's
  volume and centre by slicing the profile and integrating the clipped
  height, its area as its caps', walls' and cut faces'), all 26 cases agree
  within 2e-8 relative, every piece valid. Pieces bounded by planar faces
  only agree to 1e-15; BRepGProp's error on the pieces with elliptic faces
  (a plane across an arc wall or a circle's wall) reaches 8.8e-9
  (`disc_through_caps`).
* A plane through a vertex only, tangent to an arc wall along a ruling,
  lying in a cap or missing leaves one solid; the U profile cut across its
  prongs gives one piece below and two above.

`compare_split.py` requires every later run to reproduce these rows (on
another platform, its reviewed record) and compares the kernel's pieces with
both: the sides' totals inside the kernel's enclosures and each piece's
counts against the native ones.
