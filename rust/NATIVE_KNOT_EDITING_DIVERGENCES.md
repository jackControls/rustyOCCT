# Native knot editing observations

The bridge compares degree, periodicity, domain, every knot/multiplicity, every
control coordinate/weight and every operation's success flag. Rust must first
equal the independently solved exact fixture. Native numeric comparisons retain
the fixed `1e-10 + 2e-12*abs(expected)` budget. Removal calls native OCCT with
`1e-9` tolerance; Rust requires exact positive homogeneous preservation.

The source reference is `3d097a0328e71b826377d4814ab05ec3c3d23871`. The local
runtime is installed OCCT 7.9.3, not a build of that source revision. The first
665 sequences were captured before implementation; the checked-in capture record
retains their hashes. Two supplemental nonconstant seam cases bring the bridge
to 667 inputs. The separate fixture corpus contains 723 cases.

OCCT 7.9.3 matches 606 cases within the fixed comparison contract. The other 61
are individually pinned by version, input hash, full native bit patterns, exact
output hash and differing fields in
[`occt-knot-editing-divergences.json`](fixtures/occt-knot-editing-divergences.json).
They are not parity matches:

- 55 native insertion/removal round trips return false from removal. The exact
  coefficient system reconstructs the original positive controls and knots;
  floating insertion/inverse insertion accumulates enough error to reject them.
- Four native round trips return true but changed control values exceed the
  fixed comparison budget. Exact results equal the entire original input.
- Two aliases of a nonconstant periodic seam removal return a different cyclic
  control order. The closed geometric trace agrees, but the function of the
  original parameter is shifted. This is tracked as a parameter correspondence
  difference, not hidden inside a larger positional tolerance.

For the last case, the input is degree two, knots
`[0.5,1,2,3,4,4.5]`, all multiplicities/weights one, and poles
`[(1.5,.75,0),(1.75,1.5,.5),(1,3,2),(-1,1,1),(0,0,0)]`.
Removing either endpoint to zero produces knots `[1,2,3,4,5]`.
Exact coefficient equations and a separate Greville collocation solve recover
`[(2,1,0),(1,3,2),(-1,1,1),(0,0,0)]`. Native returns that row rotated right by one.
At parameter `u=1`, before and exact-after positions are `(1.5,2,1)`, while
native-after is `(1,.5,0)`. Native D0 reconstruction confirms those values.
Complete polynomial equality proves `native_after(u)=before(u-1)` for this case.
Rust's immutable editing contract preserves `after(u)=before(u)`.

Linux OCCT 7.6.3 matches 603 cases, with 64 separately pinned differences:
54 rejected round trips, five accepted round trips with changed control values,
three degree-25 constant-curve removals with a changed control, and the same
two seam correspondence cases. The three constant removals return
`z=4.999999999878329` for one control instead of exactly `5`; this change is
within the native removal tolerance but exceeds our fixed control comparison
budget. Rust retains the exact constant function. All 667 complete Rust output
rows are byte-identical on Linux and macOS, SHA-256
`c1f9ca57728f2cf844045d2301bc8115836a90e9776ada6a35f34347a5a061b1`.

The Linux observations came from
[CI run 35822554728](https://github.com/jackControls/rustyOCCT/actions/runs/35822554728)
at `08c97766`; that first native job correctly failed on the then-unreviewed
Linux observations. No local version review exempts another version or changed
result. Reviews never bypass Rust's independent mathematics, ordinary
regressions, or fuzzing. General
surface knot editing, degree elevation and topology integration are not covered.

```sh
python3 rust/tools/generate_knot_editing_fixtures.py --check
python3 rust/tools/compare_knot_editing.py --occt-root /path/to/occt
```
