#!/usr/bin/env python3
"""Height split and stacked fuse fixtures from an independent enumeration (M3).

split-merge-cases.txt lists scenarios: named prisms in the identity case
protocol, then operations on them (`split BODY H OP LOWER UPPER`,
`fuse A B OP OUT`, `compose I J`). split-merge-expected.tsv gives, per
scenario, each produced body's id and every entity row (id kind role ordinal
parents locator), each step's relations in canonical order or `error`, and
each composition's relations, all from split_merge_reference.py. Corpus
scenarios, built from the identity corpus with xorshift split heights, are
recorded as FNV-1a-128 digests of their newline-joined rows. No Rust result
supplies an expectation.

native_rows() writes the explicit OCCT constructions for the scenarios the
native oracle represents: one split, one fuse, or a split whose pieces are
fused back.
"""
import argparse
from dataclasses import dataclass, field, replace
from pathlib import Path

import generate_identity_fixtures as identity
from identity_reference import (encode_case, entity_text, fnv128, frame_axes, hexid, native_case,
                                number, relation_text, stored, transform_matrix)
from split_merge_reference import Rejected, compose, extruded, fuse, split

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Scenario:
    name: str
    tolerance: float
    bodies: list                     # [(name, identity Case)]
    steps: list = field(default_factory=list)


def encode(s):
    out = [f'case {s.name} {number(s.tolerance)}']
    for name, c in s.bodies:
        body = encode_case(c).splitlines()[1:-1]
        out += [f'body {name}']+body
    for step in s.steps:
        if step[0] == 'split':
            _, src, h, op, lo, hi = step
            out.append(f'split {src} {number(h)} {op} {lo} {hi}')
        elif step[0] == 'fuse':
            _, a, b, op, o = step
            out.append(f'fuse {a} {b} {op} {o}')
        else:
            out.append(f'compose {step[1]} {step[2]}')
    out.append('end')
    return '\n'.join(out)


def evaluate(s):
    """Rows: ('body', name, text), ('step', k, text), ('compose', 'i j', text)."""
    bodies = {name: extruded(c) for name, c in s.bodies}
    histories, rows, produced = [], [], []
    for k, step in enumerate(s.steps):
        try:
            if step[0] == 'split':
                _, src, h, op, lo, hi = step
                low, high, rels = split(bodies[src], h, op)
                bodies[lo], bodies[hi] = low, high
                produced += [lo, hi]
                histories.append(([bodies[src].id], [low.id, high.id], rels))
            elif step[0] == 'fuse':
                _, a, b, op, o = step
                fused, rels = fuse(bodies[a], bodies[b], op)
                bodies[o] = fused
                produced.append(o)
                histories.append(([bodies[a].id, bodies[b].id], [fused.id], rels))
            else:
                _, i, j = step
                first, second = histories[i], histories[j]
                assert first[1] == second[0], 'composition needs chained bodies'
                rels = compose(first[2], second[2])
                histories.append((first[0], second[1], rels))
                rows += [('compose', f'{i} {j}', relation_text(r)) for r in rels]
                continue
        except Rejected:
            histories.append(None)
            rows.append(('step', str(k), 'error'))
            continue
        rows += [('step', str(k), relation_text(r)) for r in rels]
    for name in produced:
        b = bodies[name]
        rows.append(('body', name, f'id {hexid(b.id)}'))
        rows += [('body', name, t) for t in sorted(entity_text(e) for e in b.entities)]
    return rows


# ------------------------------------------------------------------ scenarios

def _mid(a, b, t=0.5):
    return a+(b-a)*t


def explicit():
    out = []
    bases = [c for c in identity.explicit() if c.box is None]
    for c in bases:
        lo, hi = sorted((c.start, c.end))
        h = _mid(lo, hi, 0.375)
        h2 = _mid(h, hi, 0.5)
        # Split, fuse the pieces back, compose; split the upper piece again.
        out.append(Scenario(f'split_{c.name}', c.tolerance, [('P', c)], [
            ('split', 'P', h, 101, 'L', 'U'), ('fuse', 'L', 'U', 102, 'F'), ('compose', 0, 1),
            ('split', 'U', h2, 103, 'U0', 'U1')]))
    by = {c.name: c for c in bases}
    tol = 1e-7

    def stacked(name, base, first, second, order='ab'):
        a = replace(by[base], name='A', start=first[0], end=first[1])
        # A separate construction: its own operation id, so its own ids.
        b = replace(by[base], name='B', start=second[0], end=second[1], operation=by[base].operation+50)
        x, y = ('A', 'B') if order == 'ab' else ('B', 'A')
        return Scenario(name, tol, [('A', a), ('B', b)], [('fuse', x, y, 201, 'F')])
    out.append(stacked('fuse_polygon_5', 'polygon_5', (0.0, 2.0), (2.0, 5.0)))
    out.append(stacked('fuse_polygon_5_upper_first', 'polygon_5', (0.0, 2.0), (2.0, 5.0), 'ba'))
    out.append(stacked('fuse_decreasing', 'polygon_7', (3.0, 1.0), (1.0, -2.0)))
    out.append(stacked('fuse_circle', 'circle', (0.0, 4.0), (4.0, 4.5)))
    out.append(stacked('fuse_holes_3_labelled', 'holes_3_labelled', (0.0, 3.0), (3.0, 4.5)))
    out.append(stacked('fuse_transformed', 'transformed_copy', (-1.0, 0.5), (0.5, 2.0)))
    # Fusing the pieces in the other body order is a fuse, not an inverse.
    c = by['holes_3']
    out.append(Scenario('split_then_fuse_upper_first', tol, [('P', c)], [
        ('split', 'P', 1.25, 101, 'L', 'U'), ('fuse', 'U', 'L', 102, 'F')]))
    # Domain errors.
    p = by['polygon_6']
    for name, h in [('at_start', 0.0), ('at_end', 4.0), ('below', -1.0), ('above', 5.0),
                    ('within_tolerance', 0.5*tol), ('within_tolerance_of_end', 4.0-0.5*tol)]:
        out.append(Scenario(f'split_rejected_{name}', tol, [('P', p)], [('split', 'P', h, 101, 'L', 'U')]))
    out.append(stacked('fuse_rejected_gap', 'polygon_5', (0.0, 2.0), (2.5, 5.0)))
    out.append(stacked('fuse_rejected_overlap', 'polygon_5', (0.0, 2.0), (1.0, 5.0)))
    out.append(stacked('fuse_rejected_opposite', 'polygon_5', (0.0, 2.0), (5.0, 2.0)))
    a = replace(by['polygon_5'], name='A', start=0.0, end=2.0)
    b = replace(by['polygon_5'], name='B', start=2.0, end=4.0)
    out.append(Scenario('fuse_rejected_same_ids', tol, [('A', a), ('B', b)], [('fuse', 'A', 'B', 201, 'F')]))
    b = replace(by['polygon_6'], name='B', start=2.0, end=4.0, operation=57)
    out.append(Scenario('fuse_rejected_profile', tol, [('A', a), ('B', b)], [('fuse', 'A', 'B', 201, 'F')]))
    a = replace(by['labels_absent'], name='A', start=0.0, end=1.0)
    b = replace(by['labels_present'], name='B', start=1.0, end=2.0, operation=61)
    out.append(Scenario('fuse_rejected_labels', tol, [('A', a), ('B', b)], [('fuse', 'A', 'B', 201, 'F')]))
    a = replace(by['polygon_5'], name='A', start=0.0, end=2.0)
    b = replace(by['polygon_5'], name='B', start=2.0, end=4.0, operation=57,
                frame=(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0))
    out.append(Scenario('fuse_rejected_frame', tol, [('A', a), ('B', b)], [('fuse', 'A', 'B', 201, 'F')]))
    return out


def corpus():
    """Split and fuse back every other identity corpus prism at xorshift heights."""
    rng = identity.Xorshift(0x5b117f05e2026)
    out = []
    for c in identity.corpus()[:128]:
        t = 0.1+0.8*rng.unit()
        lo, hi = sorted((c.start, c.end))
        out.append(Scenario(f'corpus_{c.name[7:]}', c.tolerance, [('P', c)], [
            ('split', 'P', _mid(lo, hi, t), 2000, 'L', 'U'), ('fuse', 'L', 'U', 2001, 'F'),
            ('compose', 0, 1)]))
    return out


def generate():
    explicit_cases, corpus_cases = explicit(), corpus()
    scenarios = explicit_cases+corpus_cases
    names = [s.name for s in scenarios]
    assert len(names) == len(set(names))
    expected = ['# case\tsection\tname\trow (explicit) or digest (corpus)']
    for s in explicit_cases:
        expected += [f'{s.name}\t{kind}\t{key}\t{text}' for kind, key, text in evaluate(s)]
    for s in corpus_cases:
        rows = [f'{kind} {key} {text}' for kind, key, text in evaluate(s)]
        expected.append(f'{s.name}\tdigest\t{len(rows)}\t{fnv128(chr(10).join(rows).encode()).hex()}')
    return scenarios, {
        'split-merge-cases.txt': '\n'.join(encode(s) for s in scenarios)+'\n',
        'split-merge-expected.tsv': '\n'.join(expected)+'\n',
    }


# ------------------------------------------------------------------ native rows

def _apply(matrix, p, direction=False):
    return tuple(sum(matrix[i][j]*p[j] for j in range(3))+(0.0 if direction else matrix[i][3])
                 for i in range(3))


def _plane(c, h):
    """The split plane through the prism's axis point at offset h, transformed."""
    o, _, _, n = frame_axes(c.frame)
    point, normal = tuple(o[i]+n[i]*h for i in range(3)), n
    for t in c.transforms:
        m = transform_matrix(t)
        point, normal = _apply(m, point), _apply(m, normal, True)
    size = max(abs(v) for b in c.boundaries
               for v in ((b.circle[0], b.circle[1], b.circle[2]) if b.circle
                         else [x for p in stored(b, c.tolerance)[0] for x in p]))
    half = 8.0*(size+abs(c.start)+abs(c.end))
    return ' '.join(number(v) for v in (*point, *normal, half))


def _prism(name, c):
    return [f'body {name}']+native_case(c).splitlines()[1:-1]


def native_scenario(s):
    """The steps of a scenario the native oracle represents."""
    if s.steps[0][0] == 'split' and len(s.steps) >= 2 and s.steps[1][:3] == ('fuse', 'L', 'U'):
        return replace(s, steps=[s.steps[0], s.steps[1], ('compose', 0, 1)])
    return replace(s, steps=s.steps[:1])


def native_rows(scenarios):
    """(scenario, oracle input) for the scenarios the native oracle represents."""
    out = []
    for s in scenarios:
        bodies = dict(s.bodies)
        steps = s.steps
        if 'rejected' in s.name:
            continue
        if steps[0][0] == 'split' and len(steps) >= 2 and steps[1][:3] == ('fuse', 'L', 'U'):
            op = ['splitfuse P '+_plane(bodies['P'], steps[0][2])]
        elif [st[0] for st in steps] == ['split']:
            op = ['split P '+_plane(bodies['P'], steps[0][2])]
        elif [st[0] for st in steps] == ['fuse'] and len(bodies) == 2:
            op = [f'fuse {steps[0][1]} {steps[0][2]}']
        else:
            continue
        rows = [f'case {s.name}']
        for name, c in s.bodies:
            rows += _prism(name, c)
        out.append((s, '\n'.join(rows+op+['end'])))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    scenarios, files = generate()
    print(f'{len(scenarios)} scenarios, {len(native_rows(scenarios))} native')
    for name, contents in files.items():
        path = ROOT/'fixtures'/name
        if args.check:
            if path.read_text() != contents:
                parser.error(f'{name} changed; investigate before updating')
        else:
            path.write_text(contents)


if __name__ == '__main__':
    main()
