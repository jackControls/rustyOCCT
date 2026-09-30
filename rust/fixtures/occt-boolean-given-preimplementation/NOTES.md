# A stack, an S9b.1 result or one solid of several given to another Boolean (S9e.2): native observations before any kernel code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0
built headless with exception checks, captured 2026-09-30 by
`compare_given_boolean.py --capture` (`compare_boolean.py`'s probe and
protocol on S9e.2's fixtures). `capture.json` records the Rust revision
(`697d4fe4`, the evidence commit) and that the kernel's S9e.2 code did not
exist (`rust_given_boolean_exists` false:
`rust/kernel/src/solid/boolean/curved/matched.rs` absent). The kernel's
probe reports every case `unsupported`: its first Boolean evaluates (a
stack, an S9b.1 result, or two solids of which the pick point's is taken)
and the second is refused with `OutOfDomain` (a stack's or an S9b.1
result's partner with arcs by `polyhedra.rs`'s `prism_model`, or the
stack's cylinders by its `stored_model`; a result of several solids by
`curved/given.rs`), each naming S9e.2. The protocol adds a pick point: a
`then OP ID [swapped] solid X Y Z` row takes the first result's solid
holding the point (`occt_boolean_oracle.cpp`: `BRepClass3d_SolidClassifier`
at tolerance 1e-7, `TopAbs_IN` for exactly one of the result's solids, else
`failure`). `inputs.txt` holds the 36 cases of `boolean-given-cases.txt`
as explicit constructions (`identity_reference.native_chained_case`),
`oracle.cpp` is `occt_boolean_oracle.cpp` as captured and `native.txt` its
results; each prism is `BRepPrimAPI_MakePrism` of its profile face in its
own frame (`XY`, `SIDE`, `DOWN`, `TILT` and `R125`).

## Observations

* Every case is done, every result valid with no warnings, and the solid
  counts are the reference's in all 36 (`chained_boolean_reference.py`,
  its selector for the picked solids; `boss_sliced_common` two solids, the
  others one); every pick point lies inside exactly one solid of its first
  result.
* All 36 match within the 2e-8 allowance: volumes within 5.5e-9 relative,
  areas within 5.7e-9 and centres within 3.9e-10 of the case's size
  (`slanted_bored_common` and `split_bored_common`, a cylinder's piece
  under the tilted planes: BRepGProp on its ellipse-bounded faces); the
  chains in exact frames (the rollex, the boss bored, the severed box)
  within 3.3e-15. No review is needed.
* The rollex's cut (`rollex_turned_cut`, `rollex_flat_cut`) has area
  30152.95448, DRAW's `checkprops result -s 30153` in `bcut_simple/L3` to
  `L6`, both frames alike.
* The degenerate fixtures are valid results natively: the third cylinder
  internally tangent to the rollex's rim along a generatrix
  (`rollex_tangent`) and the third box whose wall touches the severed box's
  other solid along a generatrix (`severed_tangent`), all three operations
  each, as the reference measures them; the decisions refuse the first by
  S9c.1's rules (a tangency between the inputs) and the second
  conservatively (the given model's arrangement holds every solid of the
  first result).
* Unified counts: the rollex's six results and the tangent rollex common
  lose faces, edges or vertices when unified (OCCT's split faces and seams
  joined); the other solids keep their counts.
