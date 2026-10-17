#!/usr/bin/env python3
"""Fixtures for the shallow lens (HANDOFF item 5 (d), "The shallow lens
refined" in REVIEW_NOTES.md): a ball of radius 0.75 crossing a rod's wall
(radius 1, heights -3 to 3) by a depth `d` of its radius, in the level frame
and the boolean target's turned one (origin `(1, -2, 0.5)`, normal `(0, 3,
4)`, `x` hint `(1, 0, 0)`).

`shallow-lens-cases.tsv` gives each frame's stored axes and the rod's bottom
cap origin as the reference emulates `Frame3::new` and `o - 3 n` (hex bits,
which the Rust test checks), and each case's ball centre (the frame's point
at `u = 1 + 0.75 (1 - d)`, `v = 0`, height `0.25`, exact and rounded once)
with `shallow_lens_reference.py`'s 30-digit references: the lens's volume,
area and its wall's and sphere's parts, the rod's and the ball's volumes and
areas. Depths `1e-5`, `1e-6`, `3e-7` and `1.4e-7` (the last `1.05e-7` deep,
just beyond the resolution `1e-7`, which refuses shallower crossings), and
the controls `0.05` and `0.25`. No Rust result supplies a value.
"""
import argparse
from pathlib import Path
import struct

import shallow_lens_reference as ref

ROOT = Path(__file__).resolve().parents[1]
FRAMES = [
    ('level', (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
    ('turned', (1.0, -2.0, 0.5), (0.0, 3.0, 4.0), (1.0, 0.0, 0.0)),
]
DEPTHS = [1e-5, 1e-6, 3e-7, 1.4e-7, 0.05, 0.25]
RADIUS = 1.0
BALL = 0.75
HEIGHT = 0.25
KEYS = ['lens_volume', 'lens_area', 'lens_wall_area', 'lens_sphere_area',
        'rod_volume', 'rod_area', 'ball_volume', 'ball_area']


def bits(v):
    return ' '.join(struct.pack('>d', x).hex() for x in v)


def cases():
    out = []
    for name, origin, normal, hint in FRAMES:
        frame = ref.frame_axes(origin, normal, hint)
        for d in DEPTHS:
            centre = ref.world64(*frame, RADIUS+BALL*(1.0-d), 0.0, HEIGHT)
            out.append((f'{name}_{d!r}', name, frame, d, centre))
    return out


def generate():
    rows = ['# Shallow lenses: generate_shallow_lens_fixtures.py, shallow_lens_reference.py',
            '# frame\tname\torigin\tx\ty\tnormal\trod bottom origin (hex bits)']
    for name, origin, normal, hint in FRAMES:
        o, x, y, n = ref.frame_axes(origin, normal, hint)
        bottom = ref.at_height(o, n, -3.0)
        rows.append('\t'.join(['frame', name, bits(o), bits(x), bits(y), bits(n), bits(bottom)]))
    rows.append('# case\tname\tframe\tdepth (of the ball\'s radius)\tcentre (hex bits)\t'+'\t'.join(KEYS))
    for name, frame_name, frame, d, centre in cases():
        values = ref.measures(ref.Lens(frame, RADIUS, centre, BALL))
        rows.append('\t'.join(['case', name, frame_name, repr(d), bits(centre)]
                              + [ref.text(values[k]) for k in KEYS]))
    return {'shallow-lens-cases.tsv': '\n'.join(rows)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    files = generate()
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)
    print(len(cases()), 'cases')


if __name__ == '__main__':
    main()
