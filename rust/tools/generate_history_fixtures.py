#!/usr/bin/env python3
"""Complete history fixtures for every value-identity case.

For each case in identity-cases.txt, the independent enumeration in
identity_reference.py gives the construction's relations (every vertex, edge
and face Generated from its profile or meridian element with its role), each rigid
transform's relations (every entity Modified with its id) and the composed
history from the construction to the last transform. Relations are listed in
canonical order. Explicit cases are recorded row by row; the 512 corpus cases
as an FNV-1a-128 digest of their newline-joined rows. No Rust result supplies
an expectation.
"""
import argparse
from pathlib import Path

import generate_identity_fixtures as identity
from identity_reference import (extrude_entities, extrude_history, fnv128, hexid, relation_text,
                                transform_history)

ROOT = Path(__file__).resolve().parents[1]


def rows(c):
    """(step, relation text) for one case, in canonical relation order."""
    construct = [relation_text(r) for r in extrude_history(c)]
    # box_at composes its cuboid construction with a translation; ids and
    # generated relations are those of the construction.
    out = [('construct', r) for r in construct]
    for k, _ in enumerate(c.transforms):
        out += [(f'transform{k}', relation_text(r)) for r in transform_history(c)]
    if c.transforms:
        # Generated targets keep their ids through every rigid motion.
        out += [('composed', r) for r in construct]
    return out


def generate():
    cases, _ = identity.generate()
    expected = ['# case\tstep\trelation (explicit) or digest (corpus)']
    for c in cases:
        steps = rows(c)
        if c.name.startswith('corpus_'):
            # Digests cover vertex, edge and face relations; region relations
            # (added by the cell-complex migration) are listed beside them.
            regions = {hexid(e.id) for e in extrude_entities(c) if e.kind == 'region'}
            is_region = lambda r: any(i in r.split() for i in regions)
            kept = [(s, r) for s, r in steps if not is_region(r)]
            text = '\n'.join(f'{s} {r}' for s, r in kept)
            expected.append(f'{c.name}\tdigest\t{len(kept)} {fnv128(text.encode()).hex()}')
            expected += [f'{c.name}\t{s}\t{r}' for s, r in steps if is_region(r)]
        else:
            expected += [f'{c.name}\t{s}\t{r}' for s, r in steps]
    return cases, {'history-expected.tsv': '\n'.join(expected)+'\n'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases, files = generate()
    print(f'{len(cases)} cases')
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)


if __name__ == '__main__':
    main()
