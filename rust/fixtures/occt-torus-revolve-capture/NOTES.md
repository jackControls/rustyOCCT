# Torus history (S3): native MakeRevol observations before any kernel torus code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0.
`inputs.txt` is `torus-cases.txt`; `oracle.cpp` is `occt_revolve_oracle.cpp`
as captured; `native.txt` its output. `capture.json` records that no kernel
torus code and no Python enumeration of torus entities existed.

## Probe

A closed meridian is the minor disc: one closed circle edge from its point at
`angle1`. An open one is the polygon `(0, r sin a1)`, `(R + r cos a1,
r sin a1)`, the arc to `(R + r cos a2, r sin a2)` and `(0, r sin a2)`. It is
revolved with `BRepPrimAPI_MakeRevol`, a full turn when `angle` is the
binary64 `2π` and by `angle` otherwise, and queried as for the cone.

## Observations

A full turn behaves as for cones and spheres: the rim points generate the
latitude circles, the arc the torus face, and the axis, the radial edges and
the meridian face are reported deleted (the discs and the solid unreached).
The whole torus's latitude circle and meridian circle are its two seams. A
wedge is not closed: the meridian face generates the solid, its first and
last shapes are the end discs, the arc's first and last shapes the meridian
circles at both ends, and its point generates the latitude arc (a seam in
`v`).
