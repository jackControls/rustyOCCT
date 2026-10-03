#!/usr/bin/env python3
"""Source-pinned BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut and BRepAlgoAPI_Common
observations beside the S9e.4a reference: imported solids (bodies without a
construction, read from `.brep` files OCCT wrote) given to Booleans
(`generate_imported_boolean_fixtures.py`).

The protocol, the native probe (`occt_boolean_oracle.cpp`: an imported input
its `brep PATH` row, the file's one solid read by `BRepTools::Read` from
`rust/fixtures`, passed as OCCT_BOOLEAN_FIXTURES; the other inputs their
constructions) and the comparison are `compare_boolean.py`'s, run through
`compare_boolean.make_set` on S9e.4a's fixtures, capture
(`fixtures/occt-boolean-imported-preimplementation`) and reviews
(`fixtures/occt-boolean-imported-divergences.json`). The capture was taken
before the kernel's S9e.4 module (`rust/kernel/src/solid/imported.rs`)
existed, while no solid could be made from an imported topology (the probe
reads and converts every file and reports `unsupported`): the capture's
`rust_imported_boolean_exists` is false, and the comparison requires every
case `unsupported` until the module exists. With it the kernel's rows are
compared as S9c.1's, its enclosures against the reference's totals with a
slack of 1e-12 relative (`IMPORTED_SLACK`): an imported body's stored data
are OCCT's roundings of the construction the reference takes (vertices to 15
significant digits, frames normalized again by the converter), so the
kernel's model of it lies within rounding of the reference's, not on it.

`--write-bodies` writes the bodies themselves first: the oracle's `write`
blocks of `boolean-imported-bodies.txt` (OCCT's `BRepPrimAPI` makers and
`MakePrism`, `BRepTools::Write` in format version 1) into
`rust/fixtures/imported/`.
"""
import json
import math
import sys

import compare_boolean as base
import generate_imported_boolean_fixtures as fixtures

ROOT = base.ROOT
KERNEL = ROOT/'rust/kernel/src/solid/imported.rs'
FIXTURES = ROOT/'rust/fixtures'
IMPORTED_SLACK = 1e-12


def rust_imported_boolean_exists():
    """Whether the kernel's S9e.4 module exists."""
    return KERNEL.exists()


class ImportedSet(base.Set):
    def __init__(self, splines):
        assert not splines, 'S9e.4a has no spline set'
        super().__init__(False)
        self.cases = fixtures.all_cases
        self.expected = ROOT/'rust/fixtures/boolean-imported-expected.tsv'
        self.capture = ROOT/'rust/fixtures/occt-boolean-imported-preimplementation'
        self.output = ROOT/'target/boolean-imported-oracle'
        self.reviews = ROOT/'rust/fixtures/occt-boolean-imported-divergences.json'
        self.exists_key = 'rust_imported_boolean_exists'

    def exists(self):
        return rust_imported_boolean_exists()

    def before_code(self):
        return not rust_imported_boolean_exists()


def native_input():
    """The imported cases' native rows."""
    return fixtures.native_input()


def case_scale(case):
    """The case's size: its solids' reach from the origin (an imported one's
    construction's)."""
    values = [1.0]
    for c in case.constructions:
        reach = max(abs(x) for x in c.frame[:3])
        if c.sphere is not None:
            values.append(reach+c.sphere[0])
        elif c.cone is not None:
            values.append(reach+max(c.cone[0], c.cone[1])+c.cone[2])
        elif c.torus is not None:
            values.append(reach+c.torus[0]+c.torus[1])
        else:
            values += [reach+abs(c.start), reach+abs(c.end)]
            for b in c.boundaries:
                if b.circle is not None:
                    values += [reach+abs(b.circle[0])+b.circle[2], reach+abs(b.circle[1])+b.circle[2]]
                else:
                    values += [reach+abs(x) for p in b.points for x in p]
    return max(values)


def rust_differences(rust, expected, native):
    """`compare_boolean.rust_differences` with the imported bodies' slack."""
    kind = expected['expect'][0]
    if rust[0][0] == 'limit':
        return ['rust_limit']
    if rust[0][0] == 'refused':
        return ['rust_refused'] if kind == 'degenerate' else ['rust_refused_unexpectedly']
    want = expected['result']
    solids = [([float(x) for x in w[1:3]], [float(x) for x in w[3:5]],
               [[float(x) for x in w[5+2*i:7+2*i]] for i in range(3)], [int(x) for x in w[11:14]])
              for w in rust if w[0] == 'solid']
    if want is None:
        return [] if rust[0][0] == 'empty' and not solids else ['rust_solids']
    count, V, A, c = want
    if len(solids) != count:
        return ['rust_solids']
    vlo = sum(s[0][0] for s in solids)
    vhi = sum(s[0][1] for s in solids)
    alo = sum(s[1][0] for s in solids)
    ahi = sum(s[1][1] for s in solids)
    out = []
    slack = lambda x: IMPORTED_SLACK*abs(x)
    if not (vlo-slack(V) <= V <= vhi+slack(V)) or not (alo-slack(A) <= A <= ahi+slack(A)):
        out.append('rust_measure_outside_reference')
    for i in range(3):
        lo = sum(min(a*b for a in s[0] for b in s[2][i]) for s in solids)
        hi = sum(max(a*b for a in s[0] for b in s[2][i]) for s in solids)
        m = V*c[i]
        allow = IMPORTED_SLACK*max(1.0, abs(V)*math.sqrt(sum(x*x for x in c)))
        if not (lo-allow <= m <= hi+allow):
            out.append('rust_centre_outside_reference')
    if native is not None and native[0] == 'done':
        if sorted(tuple(s[3]) for s in solids) != sorted(tuple(s[4]) for s in native[3]):
            out.append('rust_entity_counts')
    return sorted(set(out))


def build(*args, **kwargs):
    """The oracle, its environment naming the fixtures directory."""
    executable, env, loaded, command = BUILD(*args, **kwargs)
    env = dict(env, OCCT_BOOLEAN_FIXTURES=str(FIXTURES))
    return executable, env, loaded, command


BUILD = base.build
base.build = build
base.make_set = ImportedSet
base.native_input = native_input
base.case_scale = case_scale
base.rust_differences = rust_differences


def write_bodies():
    """`--write-bodies`: OCCT writes every body's file."""
    import argparse
    from compare_degree_elevation import verify_sdk
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occt-root', type=base.Path, required=True)
    parser.add_argument('--sdk-manifest', type=base.Path, required=True)
    parser.add_argument('--write-bodies', action='store_true')
    args = parser.parse_args()
    output = (ROOT/'target/boolean-imported-oracle').resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = args.occt_root.resolve()
    verify_sdk(prefix, args.sdk_manifest)
    executable, env, _, _ = build(prefix, output, base.SOURCE_FILE, 'boolean-oracle',
                                  ('TKBO', 'TKPrim', 'TKShHealing'),
                                  ['TKBO', 'TKPrim', 'TKShHealing']+base.TOOLKITS)
    (FIXTURES/'imported').mkdir(exist_ok=True)
    text = (FIXTURES/fixtures.BODIES_FILE).read_text()
    record = base.run(executable, text, env)
    lines = record['stdout'].split()
    written = [w for w in record['stdout'].splitlines() if w.endswith(' written')]
    if record['exit_code'] != 0 or len(written) != len(fixtures.BODIES) or 'failure' in lines:
        raise SystemExit('the bodies were not all written: '+json.dumps(record)[:2000])
    print(f'{len(written)} bodies written')


if __name__ == '__main__':
    if '--write-bodies' in sys.argv:
        write_bodies()
    else:
        base.main()
