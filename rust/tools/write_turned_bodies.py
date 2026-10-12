#!/usr/bin/env python3
"""S9e.4b.4d's evidence: turned bodies OCCT writes (REVIEW_NOTES.md,
"S9e.4b.4d refined"), small kin of the dataset's `OCC485a` (a turned part
whose fillet and round are tori joining their neighbours smoothly, and whose
step is a cone within the resolution of flat), `OCC497c` (a pin whose step is
such a cone) and `bug25966_Revolution` (a cylinder OCCT stores as two faces on
frames apart along its axis).

Each body is one `write` block of `occt_boolean_oracle.cpp`: a profile in the
world's `xz` plane (a `plane` row, the frame's normal `-y` so that a
counter-clockwise arc in `(x, z)` is the row's flag 1) of lines and arcs, its
edge on the `z` axis, turned a whole turn about `z` by BRepPrimAPI_MakeRevol
(the oracle's `revolve` row), written by BRepTools::Write in format version 1
under `rust/fixtures/imported/`. The bodies:

* `turned_fillet`: radius 10 to height 4, a step in to a concave fillet of
  radius 0.5 (a torus tangent to the step's plane and the shaft), the shaft of
  radius 5 to height 8, a convex round of radius 2 (a torus tangent to the
  shaft and the top's plane) and the top disc of radius 3 at height 10.
* `turned_round_off`: the same but the round's centre raised by 1e-4 and its
  radius the distance to the shaft's top (2 + 2.5e-9): its tube crosses the
  shaft by 2.5e-9, as `OCC485a`'s round does by 3.76e-6, and the normals at
  the stored edge differ by 5e-5 rad.
* `turned_corner` (control): a torus meeting cones at corners (radius 7 at
  the bottom, a slanted wall to an arc of radius 2 about `(4, 4)`, a slanted
  top): no smooth join.
* `turned_pin`: radius 6 to height 5, a step in to radius 4 that rises 1e-5
  (a cone 5e-6 rad from flat), radius 4 to height 12, every tolerance set to
  1e-3 (the oracle's `tolerance` row) as a translated part stores its own: the
  step within 1 % of the resolution of a plane, as `OCC497c`'s steps (OCCT
  turns a step rising less than 1e-7 into a plane).
* `turned_step` (control): the same step rising 1e-3, ten thousand
  resolutions off flat.
* `turned_split`: radius 8 with a vertex at height 3 on its wall (two faces
  on one cylinder), a cone in to radius 5 at height 14, radius 5 to height 20.

`--check` verifies `rust/fixtures/turned-bodies.txt` against the rows here,
`--write-text` writes it, and `--write-bodies --occt-root ... --sdk-manifest
...` builds the oracle with the pinned SDK and writes the files.
"""
import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT/'rust/fixtures'
TEXT = FIXTURES/'turned-bodies.txt'
TURN = repr(2*math.pi)
PLANE = 'plane 0 0 0 0 -1 0 1 0 0'
REVOLVE = f'revolve 0 0 0 0 0 1 {TURN}'


def num(x):
    """A coordinate's shortest round-trip text."""
    return repr(float(x)) if x != int(x) else str(int(x))


def profile(points):
    """A `wire S` row of `(x, z, segment)` points, each segment `('L',)` or
    `('A', cx, cz, r, ccw)` to the next point (the last to the first)."""
    out = [f'wire S {len(points)}']
    for x, z, seg in points:
        out.append(f'{num(x)} 0 {num(z)}')
        if seg[0] == 'L':
            out.append('L')
        else:
            _, cx, cz, r, ccw = seg
            out.append(f'A {num(cx)} 0 {num(cz)} {num(r)} {1 if ccw else 0}')
    return ' '.join(out)


L = ('L',)
OFF = 1e-4
ROUND = math.sqrt(4 + OFF*OFF)

BODIES = {
    'turned_fillet': [(0, 0, L), (10, 0, L), (10, 4, L), (5.5, 4, ('A', 5.5, 4.5, 0.5, False)),
                      (5, 4.5, L), (5, 8, ('A', 3, 8, 2, True)), (3, 10, L), (0, 10, L)],
    'turned_round_off': [(0, 0, L), (10, 0, L), (10, 4, L), (5.5, 4, ('A', 5.5, 4.5, 0.5, False)),
                         (5, 4.5, L), (5, 8, ('A', 3, 8 + OFF, ROUND, True)), (3, 8 + OFF + ROUND, L),
                         (0, 8 + OFF + ROUND, L)],
    'turned_corner': [(0, 0, L), (7, 0, L), (6, 4, ('A', 4, 4, 2, True)), (4, 6, L), (0, 7, L)],
    'turned_pin': [(0, 0, L), (6, 0, L), (6, 5, L), (4, 5.00001, L), (4, 12, L), (0, 12, L)],
    'turned_step': [(0, 0, L), (6, 0, L), (6, 5, L), (4, 5.001, L), (4, 12, L), (0, 12, L)],
    'turned_split': [(0, 0, L), (8, 0, L), (8, 3, L), (8, 10, L), (5, 14, L), (5, 20, L), (0, 20, L)],
}

# The written solid's tolerances, where not OCCT's own.
TOLERANCE = {'turned_pin': 0.001}


def bodies_text():
    """The oracle's `write` blocks."""
    out = []
    for name, points in BODIES.items():
        tolerance = [f'tolerance {TOLERANCE[name]}'] if name in TOLERANCE else []
        out += [f'write {name} imported/{name}.brep', PLANE, profile(points), REVOLVE] + tolerance + ['end']
    return '\n'.join(out) + '\n'


def write_bodies(args):
    """OCCT writes every body's file."""
    sys.path.insert(0, str(ROOT/'rust/tools'))
    import compare_boolean as base
    from compare_degree_elevation import verify_sdk
    output = (ROOT/'target/boolean-turned-bodies-oracle').resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, _, _ = base.build(prefix, output, base.SOURCE_FILE, 'boolean-oracle',
                                       ('TKBO', 'TKPrim', 'TKShHealing'),
                                       ['TKBO', 'TKPrim', 'TKShHealing']+base.TOOLKITS)
    env = dict(env, OCCT_BOOLEAN_FIXTURES=str(FIXTURES))
    (FIXTURES/'imported').mkdir(exist_ok=True)
    record = base.run(executable, TEXT.read_text(), env)
    written = [w for w in record['stdout'].splitlines() if w.endswith(' written')]
    if record['exit_code'] != 0 or len(written) != len(BODIES) or 'failure' in record['stdout'].split():
        raise SystemExit('the bodies were not all written: '+json.dumps(record)[:2000])
    print(f'{len(written)} bodies written')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--write-text', action='store_true')
    parser.add_argument('--write-bodies', action='store_true')
    parser.add_argument('--occt-root', type=Path)
    parser.add_argument('--sdk-manifest', type=Path)
    args = parser.parse_args()
    if args.write_text:
        TEXT.write_text(bodies_text())
    if args.check and TEXT.read_text() != bodies_text():
        raise SystemExit(f'{TEXT.name} differs from the rows')
    if args.write_bodies:
        if args.occt_root is None or args.sdk_manifest is None:
            parser.error('--write-bodies needs --occt-root and --sdk-manifest')
        write_bodies(args)


if __name__ == '__main__':
    main()
