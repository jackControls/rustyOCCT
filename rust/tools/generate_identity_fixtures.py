#!/usr/bin/env python3
"""Value-identity fixtures from an independent encoding, hash and prism enumeration.

identity-vectors.tsv pins the version-1 derivation encoding and FNV-1a-128
digest on hand-written derivations covering every parent tag, operation kind,
entity kind and role. identity-cases.txt lists extrusion cases in a text
protocol; identity-expected.tsv gives, per case, every entity's id, role,
ordinal, parents and structural locator, computed by identity_reference.py.
The 512 prisms that mirror the invariants corpus (same xorshift stream and
structure, correctly rounded trigonometry) are recorded as an FNV-1a-128 digest
of their newline-joined sorted rows; explicit cases are recorded row by row. No Rust result
supplies an expectation.
"""
import argparse
from pathlib import Path

from brep_reference import cos_rn, sin_rn
from identity_reference import (Boundary, Case, Derivation, encode_case, entity_text,
                                extrude_entities, fnv128, hexid)

ROOT = Path(__file__).resolve().parents[1]
TAU = 6.283185307179586
X_HINT = (1.0, 0.0, 0.0)


class Xorshift:
    """The invariants test's generator (tests/invariants.rs)."""
    def __init__(self, seed):
        self.s = seed

    def unit(self):
        m = (1 << 64)-1
        self.s ^= (self.s << 13) & m
        self.s ^= self.s >> 7
        self.s ^= (self.s << 17) & m
        return (self.s >> 11)/float(1 << 53)

    def signed(self):
        return 2.0*self.unit()-1.0


def corpus():
    """Mirror of generated_prisms_preserve_geometric_invariants: 256 cases and
    their reordered, reversed-direction alternates."""
    rng = Xorshift(0xc0de5eed0cc72026)
    cases = []
    for case in range(256):
        scale = 2.0**((case % 25)-12)
        tolerance = 1e-9*scale
        count = 8+case % 13
        vertices = []
        for i in range(count):
            angle = TAU*i/count
            radius = scale*(0.8+0.6*rng.unit())
            vertices.append((radius*cos_rn(angle), radius*sin_rn(angle)))
        holes = [] if case % 2 == 0 else [Boundary(circle=(0.0, 0.0, 0.2*scale))]
        normal = (rng.signed(), rng.signed(), 0.5+rng.unit())
        origin = tuple(10.0*scale*rng.signed() for _ in range(3))
        low = -scale*(0.5+rng.unit())
        high = scale*(0.5+rng.unit())
        rng.unit()                       # cut height
        axis = (0.3+rng.unit(), rng.signed(), rng.signed())
        angle = TAU*rng.signed()
        for _ in range(32*6):            # classification samples
            rng.unit()
        frame = (*origin, *normal, *X_HINT)
        transforms = [('R', (0.0, 0.0, 0.0), axis, angle),
                      ('T', (20.0*scale, -40.0*scale, 30.0*scale))]
        cases.append(Case(f'corpus_{case}', tolerance, 1000+case, frame, low, high,
                          [Boundary(points=vertices)]+holes, transforms))
        k = 1+case % (count-1)
        reordered = (vertices[k:]+vertices[:k])[::-1]
        cases.append(Case(f'corpus_{case}_alternate', tolerance, 1000+case, frame, high, low,
                          [Boundary(points=reordered)]+holes))
    return cases


def regular(n, r=10.0, c=(0.0, 0.0)):
    return [(c[0]+r*cos_rn(TAU*k/n), c[1]+r*sin_rn(TAU*k/n)) for k in range(n)]


def square(c, h):
    return [(c[0]-h, c[1]-h), (c[0]+h, c[1]-h), (c[0]+h, c[1]+h), (c[0]-h, c[1]+h)]


def labelled(points, base):
    n = len(points)
    return Boundary(points=points, labels=(base, [base+1+j for j in range(n)],
                                           [base+101+j for j in range(n)]))


def explicit():
    tol = 1e-7
    xy = (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, *X_HINT)
    tilted = (3.0, -2.0, 5.0, 0.3, -0.4, 0.8, *X_HINT)
    cases = []
    for n in range(3, 13):
        cases.append(Case(f'polygon_{n}', tol, 7, xy, 0.0, 4.0, [Boundary(points=regular(n))]))
    cases.append(Case('circle', tol, 7, xy, 0.0, 4.0, [Boundary(circle=(1.0, 2.0, 3.0))]))
    cases.append(Case('circle_labelled', tol, 7, tilted, -1.0, 2.0,
                      [Boundary(circle=(1.0, 2.0, 3.0), labels=(9, [10], [11]))]))
    holes = [Boundary(points=square((-5.0, -5.0), 1.0)), Boundary(circle=(5.0, 5.0, 1.5)),
             Boundary(points=square((5.0, -5.0), 1.0)[::-1])]
    outer = square((0.0, 0.0), 10.0)
    for k in range(4):
        cases.append(Case(f'holes_{k}', tol, 8, tilted, 0.0, 3.0, [Boundary(points=outer)]+holes[:k]))
    lab_holes = [labelled(square((-5.0, -5.0), 1.0), 500),
                 Boundary(circle=(5.0, 5.0, 1.5), labels=(650, [651], [751])),
                 labelled(square((5.0, -5.0), 1.0)[::-1], 800)]
    cases.append(Case('holes_3_labelled', tol, 8, tilted, 0.0, 3.0, [labelled(outer, 100)]+lab_holes))
    cases.append(Case('holes_3_reversed_direction', tol, 8, tilted, 3.0, 0.0,
                      [labelled(outer, 100)]+lab_holes))
    # Clockwise input: labels follow the caller's order through the reversal.
    cw = regular(7)[::-1]
    cases.append(Case('clockwise_labelled', tol, 9, xy, 0.0, 2.0, [labelled(cw, 300)]))
    cases.append(Case('clockwise_closing_point', tol, 9, xy, 0.0, 2.0,
                      [Boundary(points=cw+[cw[0]])]))
    # Labels present, permuted and moved: ids follow labels, not positions.
    hexagon = regular(6)
    cases.append(Case('labels_absent', tol, 11, xy, 0.0, 1.0, [Boundary(points=hexagon)]))
    cases.append(Case('labels_present', tol, 11, xy, 0.0, 1.0, [labelled(hexagon, 40)]))
    perm = [3, 0, 5, 1, 4, 2]
    b = labelled(hexagon, 40)
    permuted = Boundary(points=hexagon, labels=(b.labels[0], [b.labels[1][p] for p in perm],
                                                [b.labels[2][p] for p in perm]))
    cases.append(Case('labels_permuted', tol, 11, xy, 0.0, 1.0, [permuted]))
    moved = [(x*1.25+0.5, y*0.8-0.25) for x, y in hexagon]
    cases.append(Case('labels_moved', tol, 11, xy, 0.0, 1.0, [labelled(moved, 40)]))
    inserted = hexagon[:3]+[((hexagon[2][0]+hexagon[3][0])/2, (hexagon[2][1]+hexagon[3][1])/2)]+hexagon[3:]
    lab = labelled(hexagon, 40).labels
    cases.append(Case('labels_inserted_vertex', tol, 11, xy, 0.0, 1.0, [Boundary(
        points=inserted, labels=(lab[0], lab[1][:3]+[99]+lab[1][3:], lab[2][:3]+[199]+lab[2][3:]))]))
    # Rigid copies keep every id.
    cases.append(Case('transformed_copy', tol, 12, tilted, 0.0, 2.0, [labelled(regular(5), 60)],
                      [('R', (1.0, 2.0, 3.0), (0.0, 0.6, 0.8), 1.1), ('T', (10.0, -20.0, 30.0)),
                       ('R', (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), -2.5)]))
    for sx in (1.0, -1.0):
        for sy in (1.0, -1.0):
            for sz in (1.0, -1.0):
                tag = ''.join('p' if s > 0 else 'n' for s in (sx, sy, sz))
                cases.append(Case(f'box_{tag}', tol, 0, None, 0.0, 0.0, [], [],
                                  box=((1.0, 2.0, 3.0), (4.0*sx, 5.0*sy, 6.0*sz))))
    return cases


def vectors():
    """Hand-written derivations pinning the encoding and digest."""
    e = Derivation(1, 'extrude', 'vertex', 'bottom_vertex', 0, (('label', 7),)).id()
    items = [
        ('minimal', Derivation(0, 'extrude', 'face', 'start_cap', 0, ())),
        ('label', Derivation(1, 'extrude', 'vertex', 'bottom_vertex', 0, (('label', 7),))),
        ('label_max', Derivation(2**64-1, 'extrude', 'edge', 'vertical', 0, (('label', 2**64-1),))),
        ('profile_boundary', Derivation(5, 'extrude', 'face', 'end_cap', 0,
                                        (('profile', 0, 'boundary', 0), ('profile', 3, 'boundary', 0)))),
        ('profile_segment', Derivation(5, 'extrude', 'face', 'wall', 0, (('profile', 2, 'segment', 4095),))),
        ('profile_vertex', Derivation(5, 'extrude', 'edge', 'seam', 0, (('profile', 1, 'vertex', 0),))),
        ('seam_vertex_end', Derivation(5, 'extrude', 'vertex', 'seam_vertex', 1, (('label', 3),))),
        ('entity_parent', Derivation(9, 'transform', 'edge', 'top_edge', 2**32-1, (('entity', e),))),
        ('external', Derivation(0, 'external', 'face', 'external', 17, ())),
        ('region', Derivation(8, 'extrude', 'region', 'region', 0, (('label', 1), ('profile', 1, 'boundary', 0)))),
        ('body', Derivation(42, 'extrude', 'body', 'body', 0, ())),
        ('many_parents', Derivation(3, 'extrude', 'face', 'start_cap', 0,
                                    tuple(('label', 10+k) for k in range(8)))),
        # M3: split children, cut entities and merged entities.
        ('split_child', Derivation(11, 'height_split', 'face', 'wall', 1, (('entity', e),))),
        ('cut_vertex', Derivation(11, 'height_split', 'vertex', 'cut_vertex', 0, (('entity', e),))),
        ('cut_edge', Derivation(11, 'height_split', 'edge', 'cut_edge', 1, (('entity', e),))),
        ('cut_face', Derivation(11, 'height_split', 'face', 'cut_face', 0, (('entity', e), ('entity', e[::-1])))),
        ('merged', Derivation(12, 'stacked_fuse', 'region', 'region', 0, (('entity', e), ('entity', e[::-1])))),
    ]
    rows = ['# name\tencoding hex\tid hex']
    for name, d in items:
        rows.append(f'{name}\t{d.encode().hex()}\t{hexid(d.id())}')
    return '\n'.join(rows)+'\n'


def generate():
    explicit_cases, corpus_cases = explicit(), corpus()
    cases = explicit_cases+corpus_cases
    names = [c.name for c in cases]
    assert len(names) == len(set(names))
    expected = ['# case\tid kind role ordinal parents locator (explicit) or digest (corpus)']
    for c in explicit_cases:
        for row in sorted(entity_text(e) for e in extrude_entities(c)):
            expected.append(f'{c.name}\t{row}')
    # Corpus digests cover vertices, edges and faces; regions (added by the
    # cell-complex migration) are listed beside them.
    for c in corpus_cases:
        entities = extrude_entities(c)
        rows = sorted(entity_text(e) for e in entities if e.kind != 'region')
        digest = fnv128('\n'.join(rows).encode()).hex()
        expected.append(f'{c.name}\tdigest {len(rows)} {digest}')
        expected += [f'{c.name}\t{entity_text(e)}' for e in entities if e.kind == 'region']
    return cases, {
        'identity-vectors.tsv': vectors(),
        'identity-cases.txt': '\n'.join(encode_case(c) for c in cases)+'\n',
        'identity-expected.tsv': '\n'.join(expected)+'\n',
    }


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
