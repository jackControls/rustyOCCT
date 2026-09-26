# Sphere history (S3): native MakeRevol observations before any kernel sphere code

Source reference `3d097a0328e71b826377d4814ab05ec3c3d23871`, OCCT 8.1.0.
`inputs.txt` is `sphere-cases.txt`; `oracle.cpp` is
`occt_revolve_oracle.cpp` as captured; `native.txt` its output.
`capture.json` records that no kernel sphere code and no Python enumeration
of sphere entities existed (R7 of `REVIEW_NOTES.md`, unlike the cone's R11).

## Probe

Each sphere's meridian is the face bounded by `(0, R sin a1)`,
`(R cos a1, R sin a1)`, the arc of radius `R` about the origin through the
frame's `x` to `(R cos a2, R sin a2)`, and `(0, R sin a2)`, in (radius,
height); at a pole the rim point is the axis point and the axis point is
left out. It is revolved a full turn with `BRepPrimAPI_MakeRevol`, and
`BRepTools_History` and `FirstShape`/`LastShape` are queried for every
meridian subshape, as for the cone (`../occt-revolve-history-capture`).

## Observations

Every revolved solid is BRepCheck-valid with the counts `MakeSphere` gives.
The rim points generate the latitude circles or, at a pole, a degenerated
edge; the arc generates the spherical face. As for the cone, the axis, the
radial edges and the meridian face are reported deleted, and the discs and
the solid are reached by no query. A full sphere's poles generate only
degenerated edges; their vertices are used by nothing else.
