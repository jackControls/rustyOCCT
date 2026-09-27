#!/usr/bin/env python3
"""The data/occ .brep corpus as brep_io_reference.py reads it (T2).

brep-io-expected.tsv has, per file, a `file NAME` row with the unsupported
geometry records by name and count, then a `solid NAME RECORD VERDICT V E W F
SH SO` row per solid reached from the root: VERDICT is `representable` when
the structure is within the kernel's model and `unsupported` otherwise, and
the counts are OCCT's distinct subshapes of the original seamed solid. Since
S6 a `free NAME RECORD KIND VERDICT V E W F SH SO` row follows for every
shell, face, wire, edge or vertex reached from the root outside a solid
(KIND `Sh`, `Fa`, `Wi`, `Ed` or `Ve`). No Rust result supplies an
expectation.
"""
import argparse
from pathlib import Path

from brep_io_reference import files, free_shapes, summary

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT.parent/'data'/'occ'


def rows():
    out = []
    for path in files(CORPUS):
        unsupported, solids = summary(path.read_text())
        names = ','.join(f'{k}={unsupported[k]}' for k in sorted(unsupported)) or '-'
        out.append(f'file\t{path.name}\t{names}')
        for record, ok, counts in solids:
            verdict = 'representable' if ok else 'unsupported'
            out.append('\t'.join(['solid', path.name, str(record), verdict]+[str(c) for c in counts]))
        for record, kind, ok, counts in free_shapes(path.read_text()):
            verdict = 'representable' if ok else 'unsupported'
            out.append('\t'.join(['free', path.name, str(record), kind, verdict]+[str(c) for c in counts]))
    return '\n'.join(out)+'\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    target = ROOT/'fixtures'/'brep-io-expected.tsv'
    text = rows()
    if args.check:
        if target.read_text() != text:
            raise SystemExit(f'{target.name} is stale')
        print(f'{target.name} matches')
    else:
        target.write_text(text)


if __name__ == '__main__':
    main()
