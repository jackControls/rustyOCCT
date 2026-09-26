"""Independent reference for attributes through operations (M4 of
IDENTITY_AND_HISTORY.md, contract 4), built on the M3 reference's relations
without Rust.

A policy is (on_modify, on_transform, on_split, on_merge). Every input
attribute yields one outcome (key, from, to, result), derived from the one
relation its entity is a source of:

    unchanged            kept, to the same id
    modified (transform) on_transform: keep -> kept; drop -> dropped;
                         recompute -> recomputed with the callback's value
    modified (otherwise) on_modify: keep -> kept; drop -> dropped
    split                on_split: copy -> copied to every child; drop
    merged               on_merge: keep_if_equal -> kept when every parent
                         carries the key with the same value, else conflict
                         (a missing value counts as unequal); drop
    deleted              dropped

Kept, copied and recomputed values land on their targets; dropped and
conflicting ones land nowhere. Outcomes sort by key, then input id.
The fixtures' recompute callback appends the byte 0x2a to the old value.
"""
from itertools import product

MODIFY = ('keep', 'drop')
TRANSFORM = ('keep', 'drop', 'recompute')
SPLIT = ('copy', 'drop')
MERGE = ('keep_if_equal', 'drop')
POLICIES = list(product(MODIFY, TRANSFORM, SPLIT, MERGE))


def recompute(value):
    return value+b'\x2a'


def outcomes(relations, before, policies, transform):
    """(outcomes, after) for one step: `before` maps id -> {key: value};
    `after` holds the attributes that land on output ids."""
    source = {}
    for r in relations:
        kind = r[0]
        if kind == 'generated':
            continue
        sources = r[1] if kind == 'merged' else (r[1],)
        for s in sources:
            source[s] = r
    out, after = [], {}
    for fid in sorted(before):
        for key in sorted(before[fid]):
            value = before[fid][key]
            modify, on_transform, split, merge = policies[key]
            r = source[fid]
            kind = r[0]
            if kind == 'unchanged':
                to, result = [r[2]], 'kept'
            elif kind == 'modified' and transform:
                if on_transform == 'keep':
                    to, result = [r[2]], 'kept'
                elif on_transform == 'drop':
                    to, result = [], 'dropped'
                else:
                    to, result, value = [r[2]], 'recomputed', recompute(value)
            elif kind == 'modified':
                to, result = ([r[2]], 'kept') if modify == 'keep' else ([], 'dropped')
            elif kind == 'split':
                to, result = (list(r[2]), 'copied') if split == 'copy' else ([], 'dropped')
            elif kind == 'merged':
                if merge == 'drop':
                    to, result = [], 'dropped'
                elif all(before.get(p, {}).get(key) == value for p in r[1]):
                    to, result = [r[2]], 'kept'
                else:
                    to, result = [], 'conflict'
            else:
                to, result = [], 'dropped'
            for t in to:
                after.setdefault(t, {})[key] = value
            out.append((key, fid, tuple(to), result))
    # One outcome per (key, input id): that pair orders them.
    out.sort(key=lambda o: (o[0], o[1]))
    return out, after


def outcome_text(o):
    key, fid, to, result = o
    return f'outcome {key} {fid.hex()} {",".join(t.hex() for t in to) or "-"} {result}'
