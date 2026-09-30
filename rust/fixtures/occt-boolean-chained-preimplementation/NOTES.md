# A Boolean's result given to another Boolean (S9e.1): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_chained_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9e.1's fixtures). `capture.json` records the Rust revision
(`ff39281e`, the evidence commit) and that the kernel's S9e.1 code did not
exist (`rust_chained_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/given.rs` absent). The kernel's probe
reports every case `unsupported`: its first Boolean (S9c.1's, one solid in
every case) evaluates and the second is refused with `OutOfDomain("a
Boolean of a solid with curved faces or edges in any position (S9c)")`
(S9b.2's stored model, reached for the first result and, where the third
prism has an arc, for the third prism too). The protocol is new: a `then
OP` row (`then OP swapped`) and a third prism's rows chain a second Boolean
on the first result's one solid (`occt_boolean_oracle.cpp`: the first
result's solids explored, exactly one required, the second Boolean run on
it and the third prism, swapped: the third prism the argument); the rows
the oracle prints are the second Boolean's. `inputs.txt` holds the 30 cases
of `boolean-chained-cases.txt` as explicit constructions
(`identity_reference.native_chained_case`), `oracle.cpp` is
`occt_boolean_oracle.cpp` as captured and `native.txt` its results; each
prism is `BRepPrimAPI_MakePrism` of its profile face in its own frame
(`XY`, `SIDE` and `TILT`).

## Observations

* Every case is done, every result valid with no warnings, and the solid
  counts are the reference's in all 30 (`chained_boolean_reference.py`;
  `hole_tool_cut` two solids, `hole_capped_common` none, the others one).
* All 30 match within the 2e-8 allowance: volumes within 5.7e-10
  relative and areas within 8.0e-10 (`pin_step_common`, the pin's
  ellipse-bounded faces), centres within 5.9e-10 of the case's size
  (`pin_slab_common`); the chains in exact frames (`groove_drilled`,
  `quarter_bored`, `groove_tangent`, `groove_edge`) within 2.9e-16. No
  review is needed.
* The degenerate fixtures are valid results natively: the third cylinder
  tangent to the groove's wall along a generatrix (`groove_tangent`) and the
  third box whose top edge lies on the groove's wall (`groove_edge`), all
  three operations each, as the reference measures them; the decisions refuse
  them by S9c.1's rules (a tangency between the inputs, an edge of one input
  on the other's face).
* Unified counts: seven solids lose faces, edges or vertices when unified
  (OCCT's split faces and seams joined): the `hole_halved`, `hole_tool`
  and `pin_step` commons, one of `hole_tool_cut`'s two solids, the
  `hole_capped` cut and the `groove_tangent` fuse and common; the others
  keep their counts.
