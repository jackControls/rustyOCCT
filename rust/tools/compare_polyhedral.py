#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9b reference: Booleans of polyhedral prisms in any
relative position (`generate_polyhedral_fixtures.py`,
`polyhedral_reference.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`) and the
comparison are `compare_boolean.py`'s, run on S9b's fixtures, capture
(`fixtures/occt-boolean-polyhedra-preimplementation`: taken before the
kernel's S9b module `rust/kernel/src/solid/boolean/polyhedra.rs` exists) and
reviews (`fixtures/occt-boolean-polyhedra-divergences.json`). Until S9b the
kernel refuses frames with different axes, so every case is listed under
`rust_unsupported`.
"""
import compare_boolean as base
import generate_polyhedral_fixtures as fixtures

ROOT = base.ROOT
base.fixtures = fixtures
base.CAPTURE = ROOT/'rust/fixtures/occt-boolean-polyhedra-preimplementation'
base.REVIEWS = ROOT/'rust/fixtures/occt-boolean-polyhedra-divergences.json'
base.KERNEL_FILE = ROOT/'rust/kernel/src/solid/boolean/polyhedra.rs'


def expected_rows():
    """{case: {'expect': (kind, step), 'result': (N, volume, area, centre) or None, 'slabs': []}}."""
    out = {}
    for line in (ROOT/'rust/fixtures/boolean-polyhedra-expected.tsv').read_text().splitlines()[1:]:
        name, row = line.split('\t')
        w = row.split()
        e = out.setdefault(name, {'result': None, 'slabs': []})
        if w[0] == 'expect':
            e['expect'] = (w[1], w[2])
        elif w[0] == 'result':
            e['result'] = (int(w[1]), float(w[2]), float(w[3]), [float(x) for x in w[4:7]])
    return out


base.expected_rows = expected_rows

if __name__ == '__main__':
    base.main()
