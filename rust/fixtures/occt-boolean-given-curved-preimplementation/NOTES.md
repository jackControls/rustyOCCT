# Given results of spheres, cones and tori, deeper chains and given results against a sphere, cone or torus (S9e.3a): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-10-02 by
`compare_given_curved_boolean.py --capture` (`compare_boolean.py`'s probe
and protocol on S9e.3a's fixtures). `capture.json` records the Rust
revision (`81c67f77`, the evidence commit) and that the kernel's S9e.3a
code did not exist (`rust_given_curved_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/chain.rs` absent). The kernel's probe
reports every case `unsupported`: its first Boolean evaluates (and, in a
deeper chain, its second) and the Boolean given the result is refused with
`OutOfDomain` by `curved/given.rs` (a result of solids other than prisms)
or `curved/mod.rs` (a given result against a sphere, cone or torus), each
naming S9e.3. The protocol chains further Booleans: a case may hold several
`then` rows, each followed by its solid's rows, the previous Boolean's one
solid the next one's argument (`occt_boolean_oracle.cpp`'s stages; the
output is the last Boolean's). `inputs.txt` holds the 48 cases of
`boolean-given-curved-cases.txt` as explicit constructions
(`identity_reference.native_chained_case` with its further stages),
`oracle.cpp` is `occt_boolean_oracle.cpp` as captured and `native.txt` its
results; each prism is `BRepPrimAPI_MakePrism` of its profile face in its
own frame (`XY`, `SIDE`, `TILT` and `R125`), each sphere, cone and torus its
`BRepPrimAPI` primitive.

## Observations

* Every case is done, every result valid, and the solid counts are the
  reference's in all 48 (`chained_curved_boolean_reference.py`, its counts
  the declared ones by rays at two resolutions; `dome_tilt_cut`,
  `dome_drill_cut`, `groove_tilt_cut`, `holes_tilt_cut` and the touching
  `dome_touch_fuse` two solids, `dome_touch_common` empty, the others one).
  The six results with the groove's torus report warnings
  (`HasWarnings`); they are valid.
* 33 match within the 2e-8 allowance (volumes within 1.6e-8 relative,
  areas within 1.0e-8 and centres within 2.1e-9 of the case's size, the
  groove's commons under the torus; the coaxial chains in exact frames
  within 1e-15). 15 are reviewed
  (`occt-boolean-given-curved-divergences.json`): BRepGProp's default
  integration on faces bounded by the Boolean's approximated sections of
  two curved faces (a drill's or a peg's meeting with a sphere, a drill's
  with the groove's torus, the hole's cylinder's with a sphere, a cone and
  a torus) misses by up to 6.4e-6 in volume (`holed_torus_common`); the
  same native results measured with an adaptive accuracy of 1e-10 (a
  diagnostic build, not the capture's probe) are within 2.4e-9 of the
  reference, unchanged at 1e-12.
* DRAW's `bcut_simple/G9` and `H3` (`g9_cut`): the native cut's area is
  727.4813665, DRAW's `checkprops result -s 727.481`; with the
  `pcylinder` moved clear of the frustum's top circle (`g9_clear`) its
  measures are the same (the rod inside the frustum's top disc either
  way). The decisions refuse `G9` itself as a tangency (the rod's circle
  touches the top circle at `(6, 0, 4)`), and the reference measures it as
  the regularized set.
* The degenerate fixtures are valid results natively: `G9`'s three
  operations, and a box whose bottom face touches the dome's top point
  (`dome_touch`: the fuse two solids touching at the point, the common
  empty).
* Unified counts: 18 results lose faces, edges or vertices when unified
  (OCCT's split faces and seams joined, the G9 body's split bottom disc
  and wall among them); the others keep their counts.
