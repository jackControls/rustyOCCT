#!/usr/bin/env python3
"""Regenerate complete exact degree-elevation grids, independently of Rust."""
import argparse
from pathlib import Path
from degree_elevation_reference import elevate, elevate_surface
import knot_editing_reference as curve_protocol
import surface_knot_reference as surface_protocol

ROOT = Path(__file__).resolve().parents[2]


def generate_curves():
    rows = []
    for line in (ROOT/'rust/fixtures/degree-elevation-curves.txt').read_text().splitlines():
        name, curve, operations = curve_protocol.parse(line)
        flags = []
        for op, _u, degree in operations:
            assert op == 'D'
            curve = elevate(curve, degree)
            flags.append(True)
        rows.append(line+' | '+curve_protocol.encode(name, flags, curve))
    return '\n'.join(rows)+'\n'


def generate_surfaces():
    rows = []
    for line in (ROOT/'rust/fixtures/degree-elevation-surfaces.txt').read_text().splitlines():
        name, surface, operations = surface_protocol.parse(line)
        flags = []
        for op, axis, u, degree in operations:
            assert op == 'D' and axis == 2 and u.denominator == 1
            result = elevate_surface(surface, int(u), degree)
            assert result == elevate_surface(surface, int(u), degree, order=(1, 0))
            assert tuple(map(surface_protocol.domain, result[0])) == tuple(map(surface_protocol.domain, surface[0]))
            surface = result
            flags.append(True)
        rows.append(line+' | '+surface_protocol.encode(name, flags, surface, compact=True))
    return '\n'.join(rows)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for kind, generate in [('curves', generate_curves), ('surfaces', generate_surfaces)]:
        data = generate()
        path = ROOT/f'rust/fixtures/degree-elevation-{kind}.tsv'
        if args.check:
            if not path.exists() or path.read_text() != data:
                raise SystemExit(f'{path}: fixture differs from independent coefficient equations')
        else:
            path.write_text(data)
        print(f'{len(data.splitlines())} complete {kind} coefficient grids verified')


if __name__ == '__main__':
    main()
