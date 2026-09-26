#!/usr/bin/env python3
"""Attribute fixtures from an independent enumeration (M4).

attribute-cases.txt lists every explicit split/fuse scenario of M3 whose
steps succeed, each initial body first moved by a translation (a
`transform BODY OP DX DY DZ OUT` step) so that transform policies apply, with
48 attribute keys: `policy KEY modify transform split merge` for each of the
24 policy combinations, twice. Odd keys give each entity a value from its
locator alone (split pieces and their fused whole agree); even keys add the
body's name (bodies built independently disagree); a value is the first 8
bytes of that text's FNV-1a-128 digest. `attribute BODY KEY HEX LOCATOR` rows
put key k on an entity when (7i + 3k + b) mod 4 = 0 for the entity's index i
in its body and the body's index b, so about a quarter of the keys sit on
each entity and some merged parents lack a key. attribute-expected.tsv gives each step's relations and outcomes
and every produced body's attributes, from attribute_reference.py; two
scenarios row by row, the rest as FNV-1a-128 digests. No Rust result
supplies an expectation.
"""
import argparse
from dataclasses import replace
from pathlib import Path

from attribute_reference import POLICIES, outcome_text, outcomes
from generate_split_merge_fixtures import explicit as split_merge_scenarios
from identity_reference import encode_case, fnv128, number, relation_sort_key, relation_text
from split_merge_reference import Body, Rejected, extruded, fuse, split

ROOT = Path(__file__).resolve().parents[1]
SHIFT = (1.0, -2.0, 0.5)
KEYS = {1+2*k+scheme: POLICIES[k] for k in range(len(POLICIES)) for scheme in (0, 1)}
EXPLICIT = ('split_polygon_4', 'fuse_polygon_5_upper_first')


def scenarios():
    out = []
    for s in split_merge_scenarios():
        if 'rejected' in s.name:
            continue
        names = [n for n, _ in s.bodies]
        steps = [('transform', n, 90+i, SHIFT, n+'t') for i, n in enumerate(names)]
        moved = {n: n+'t' for n in names}
        for st in s.steps:
            if st[0] == 'split':
                steps.append(('split', moved.get(st[1], st[1]), *st[2:]))
            elif st[0] == 'fuse':
                steps.append(('fuse', moved.get(st[1], st[1]), moved.get(st[2], st[2]), *st[3:]))
        out.append(replace(s, name='attributes_'+s.name[:], steps=steps))
    return out


def initial_attributes(s):
    """{body: {id: {key: value}}} by the documented rule."""
    out = {}
    for b, (name, c) in enumerate(s.bodies):
        body = extruded(c)
        attrs = {}
        for i, e in enumerate(body.entities):
            loc = ' '.join(map(str, e.locator))
            for key in KEYS:
                if (7*i+3*key+b) % 4 != 0:
                    continue
                text = f'{key}:{loc}' if key % 2 == 1 else f'{key}:{loc}:{name}'
                attrs.setdefault(e.id, {})[key] = fnv128(text.encode())[:8]
        out[name] = (body, attrs)
    return out


def encode(s):
    rows = [f'case {s.name} {number(s.tolerance)}']
    for name, c in s.bodies:
        rows += [f'body {name}']+encode_case(c).splitlines()[1:-1]
    for key, (m, t, sp, mg) in KEYS.items():
        rows.append(f'policy {key} {m} {t} {sp} {mg}')
    for name, (body, attrs) in initial_attributes(s).items():
        loc = {e.id: ' '.join(map(str, e.locator)) for e in body.entities}
        for fid in sorted(attrs):
            for key in sorted(attrs[fid]):
                rows.append(f'attribute {name} {key} {attrs[fid][key].hex()} {loc[fid]}')
    for st in s.steps:
        if st[0] == 'transform':
            _, b, op, v, o = st
            rows.append(f'transform {b} {op} '+' '.join(number(x) for x in v)+f' {o}')
        elif st[0] == 'split':
            _, b, h, op, lo, hi = st
            rows.append(f'split {b} {number(h)} {op} {lo} {hi}')
        else:
            _, a, b, op, o = st
            rows.append(f'fuse {a} {b} {op} {o}')
    return '\n'.join(rows+['end'])


def evaluate(s):
    """(section, key, text) rows: each step's relations then outcomes, then
    every produced body's attributes."""
    state = initial_attributes(s)
    bodies = {n: b for n, (b, _) in state.items()}
    attrs = {n: a for n, (_, a) in state.items()}
    rows, produced = [], []
    for k, st in enumerate(s.steps):
        if st[0] == 'transform':
            _, b, _, v, o = st
            body = bodies[b]
            moved = Body(replace(body.case, transforms=body.case.transforms+[('T', v)]), body.id, body.entities)
            rels = sorted((('modified', e.id, e.id) for e in body.entities), key=relation_sort_key)
            outs, after = outcomes(rels, attrs[b], KEYS, True)
            bodies[o], attrs[o] = moved, after
            produced.append(o)
        elif st[0] == 'split':
            _, b, h, op, lo, hi = st
            lower, upper, rels = split(bodies[b], h, op)
            outs, after = outcomes(rels, attrs[b], KEYS, False)
            for name, piece in ((lo, lower), (hi, upper)):
                ids = {e.id for e in piece.entities}
                bodies[name], attrs[name] = piece, {i: v for i, v in after.items() if i in ids}
            produced += [lo, hi]
        else:
            _, a, b, op, o = st
            fused, rels = fuse(bodies[a], bodies[b], op)
            outs, after = outcomes(rels, {**attrs[a], **attrs[b]}, KEYS, False)
            bodies[o], attrs[o] = fused, after
            produced.append(o)
        rows += [('step', str(k), relation_text(r)) for r in rels]
        rows += [('step', str(k), outcome_text(x)) for x in outs]
    for name in produced:
        lines = sorted(f'{fid.hex()} {key} {value.hex()}'
                       for fid, keyed in attrs[name].items() for key, value in keyed.items())
        rows += [('attr', name, line) for line in lines]
    return rows


def generate():
    all_scenarios = scenarios()
    expected = ['# case\tsection\tname\trow (explicit) or digest']
    for s in all_scenarios:
        rows = evaluate(s)
        if s.name[len('attributes_'):] in EXPLICIT:
            expected += [f'{s.name}\t{a}\t{b}\t{c}' for a, b, c in rows]
        else:
            text = '\n'.join(f'{a} {b} {c}' for a, b, c in rows)
            expected.append(f'{s.name}\tdigest\t{len(rows)}\t{fnv128(text.encode()).hex()}')
    return all_scenarios, {
        'attribute-cases.txt': '\n'.join(encode(s) for s in all_scenarios)+'\n',
        'attribute-expected.tsv': '\n'.join(expected)+'\n',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    scenarios_, files = generate()
    print(f'{len(scenarios_)} scenarios, {len(KEYS)} keys')
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)


if __name__ == '__main__':
    main()
