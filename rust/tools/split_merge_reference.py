"""Independent reference for the height split and the stacked fuse (M3 of
IDENTITY_AND_HISTORY.md), built on identity_reference.py without Rust or OCCT.

A body is a prism: its case (profile, frame, offsets, transforms), its body id
and every entity with its derivation and structural locator. The locator side
names ('start', 'end') refer to the body's own extrusion offsets, which keep
the direction of the prism they came from.

Height split at offset h, strictly inside the prism, with operation id `op`:
    pieces in axial order (increasing offset along the frame normal):
    ordinal 0 is the lower piece, 1 the upper.
    * entities on the lower side (cap, cap edges, cap vertices) stay
      `Unchanged` in the lower piece, those on the upper side in the upper;
    * walls, vertical edges and the solid region `Split` into one child per
      piece: Derivation(op, height_split, kind, same role, piece ordinal,
      [Entity(parent)]);
    * each piece's cut side is `Generated`: the cut face (role cut_face) from
      every wall in profile order, each cut edge (cut_edge) from the wall of
      its segment, each cut vertex (cut_vertex) from the vertical edge at its
      vertex, all with the piece ordinal;
    * the input body is replaced by two bodies,
      Derivation(op, height_split, body, body, piece ordinal, [Entity(input)]).

Stacked fuse of bodies a and b, which must have the same tolerance, profile
(points and labels), frame and transforms, the same direction and adjacent
offsets (a.end == b.start or b.end == a.start). Parents are in axial order
(the lower body first), so the result does not depend on the call order:
    * walls, vertical edges and regions pair by locator and are `Merged`
      into Derivation(op, stacked_fuse, kind, same role, 0,
      [Entity(from lower), Entity(from upper)]);
    * the shared cap sides of both bodies (caps, cap edges, cap vertices)
      are `Deleted`;
    * the outer sides stay `Unchanged`;
    * the body is Derivation(op, stacked_fuse, body, body, 0, [lower, upper]);
    * bodies that share an id (the same construction twice) are rejected.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F

from identity_reference import (Case, Derivation, Entity, body_id, extrude_entities, relation_sort_key,
                                stored)


class Rejected(Exception):
    """The operation is outside its domain; the Rust call must return an error."""


@dataclass
class Body:
    case: Case            # the prism's construction (offsets are the body's own)
    id: bytes
    entities: list        # Entity, in slot order


def extruded(c):
    return Body(c, body_id(c.operation), extrude_entities(c))


def _low_high(c):
    return (c.start, c.end) if c.start < c.end else (c.end, c.start)


def _side_of(c, locator):
    """'low' or 'high' for an entity on one end of the prism, else None."""
    side = locator[-1] if locator[0] != 'region' else None
    if side not in ('start', 'end'):
        return None
    low_name = 'start' if c.start < c.end else 'end'
    return 'low' if side == low_name else 'high'


def _swept(e):
    """Walls, vertical edges and the region: the entities that span the prism."""
    return e.locator[0] == 'region' or e.locator[-1] == 'both'


def _walls(body):
    return [e for e in body.entities if e.kind == 'face' and e.locator[-1] == 'both']


def split(body, h, op):
    """(lower, upper, relations) of Solid::split_at_height."""
    c = body.case
    low, high = _low_high(c)
    if not (F(low) < F(h) < F(high)):
        raise Rejected('split height outside the prism')
    for a, b in ((low, h), (h, high)):
        # Rust compares the rounded height with the tolerance, as extrusion does.
        if not (b-a > c.tolerance):
            raise Rejected('degenerate piece')
    increasing = c.start < c.end
    offsets = [(c.start, h) if increasing else (h, c.end),
               (h, c.end) if increasing else (c.start, h)]
    walls = [('entity', w.id) for w in _walls(body)]
    by_locator = {e.locator: e for e in body.entities}
    pieces, relations = [], []
    children = {e.id: [] for e in body.entities if _swept(e)}
    for k, keep in enumerate(('low', 'high')):
        start, end = offsets[k]
        case = replace(c, start=start, end=end)
        ents = []
        for e in body.entities:
            side = _side_of(c, e.locator)
            d = e.derivation
            if _swept(e):
                child = Entity(e.kind, Derivation(op, 'height_split', e.kind, d.role, k, (('entity', e.id),)),
                               e.locator)
                ents.append(child)
                children[e.id].append(child.id)
            elif side == keep:
                ents.append(e)
            else:
                # The cut side of this piece.
                if e.kind == 'face':
                    role, parents = 'cut_face', tuple(walls)
                elif e.kind == 'edge':
                    b, _, j, _ = e.locator
                    role, parents = 'cut_edge', (('entity', by_locator[(b, 'segment', j, 'both')].id),)
                else:
                    b, _, j, _ = e.locator
                    role, parents = 'cut_vertex', (('entity', by_locator[(b, 'vertex', j, 'both')].id),)
                cut = Entity(e.kind, Derivation(op, 'height_split', e.kind, role, k, parents), e.locator)
                ents.append(cut)
                relations.append(('generated', parents, cut.id, role))
        pieces.append(Body(case, body_id(op, 'height_split', [('entity', body.id)], k), ents))
    for e in body.entities:
        if _swept(e):
            relations.append(('split', e.id, tuple(children[e.id])))
        else:
            relations.append(('unchanged', e.id, e.id))
    _unique(pieces)
    return pieces[0], pieces[1], sorted(relations, key=relation_sort_key)


def _same_profile(a, b):
    tol = a.tolerance
    if len(a.boundaries) != len(b.boundaries):
        return False
    for p, q in zip(a.boundaries, b.boundaries):
        if (p.circle is None) != (q.circle is None) or stored(p, tol) != stored(q, tol):
            return False
        if p.circle is not None and p.circle != q.circle:
            return False
    return True


def fuse(a, b, op):
    """(fused, relations) of Solid::fuse_stacked(a, b)."""
    ca, cb = a.case, b.case
    if (ca.tolerance != cb.tolerance or ca.frame != cb.frame or ca.transforms != cb.transforms
            or ca.box is not None or cb.box is not None or not _same_profile(ca, cb)):
        raise Rejected('different prisms')
    if (ca.start < ca.end) != (cb.start < cb.end):
        raise Rejected('opposite directions')
    # I4: ids name one entity in a lineage; bodies that share one (the same
    # construction twice) cannot be combined.
    if {e.id for e in a.entities} & {e.id for e in b.entities} or a.id == b.id:
        raise Rejected('id collision')
    if ca.end == cb.start:
        start, end = ca.start, cb.end
    elif cb.end == ca.start:
        start, end = cb.start, ca.end
    else:
        raise Rejected('not adjacent')
    case = replace(ca, start=start, end=end)
    lower, upper = (a, b) if _low_high(ca)[0] < _low_high(cb)[0] else (b, a)
    in_upper = {e.locator: e for e in upper.entities}
    ents, relations = [], []
    for e in lower.entities:
        if _swept(e):
            other = in_upper[e.locator]
            d = e.derivation
            merged = Entity(e.kind, Derivation(op, 'stacked_fuse', e.kind, d.role, 0,
                                               (('entity', e.id), ('entity', other.id))), e.locator)
            ents.append(merged)
            relations.append(('merged', (e.id, other.id), merged.id))
    for body, keep in ((lower, 'low'), (upper, 'high')):
        for e in body.entities:
            if _swept(e):
                continue
            if _side_of(body.case, e.locator) == keep:
                ents.append(e)
                relations.append(('unchanged', e.id, e.id))
            else:
                relations.append(('deleted', e.id))
    fused = Body(case, body_id(op, 'stacked_fuse', [('entity', lower.id), ('entity', upper.id)]), ents)
    _unique([fused])
    return fused, sorted(relations, key=relation_sort_key)


def compose(first, second):
    """The relations of `first` followed by `second` (H5), for chains where
    every entity flows one-to-one, one-to-many or many-to-one."""
    flow1, flow2, gen1, gen2 = {}, {}, [], []
    for flows, gens, rels in ((flow1, gen1, first), (flow2, gen2, second)):
        for r in rels:
            if r[0] == 'generated':
                gens.append(r)
            elif r[0] == 'deleted':
                flows[r[1]] = []
            elif r[0] in ('unchanged', 'modified'):
                flows[r[1]] = [r[2]]
            elif r[0] == 'split':
                flows[r[1]] = list(r[2])
            else:
                for f in r[1]:
                    flows[f] = [r[2]]
    finals = {s: [x for t in ts for x in flow2[t] if True] for s, ts in flow1.items()}
    for s in finals:
        out = []
        for x in finals[s]:
            if x not in out:
                out.append(x)
        finals[s] = out
    contributors = {}
    for s, out in finals.items():
        for x in out:
            contributors.setdefault(x, []).append(s)
    unchanged1 = {r[1] for r in first if r[0] == 'unchanged'}
    unchanged2 = {r[1] for r in second if r[0] == 'unchanged'}
    rels, merged = [], set()
    for s, out in finals.items():
        if not out:
            rels.append(('deleted', s))
        elif len(out) == 1:
            x = out[0]
            if len(contributors[x]) == 1:
                if s in unchanged1 and s in unchanged2 and x == s:
                    rels.append(('unchanged', s, s))
                else:
                    rels.append(('modified', s, x))
            elif x not in merged:
                merged.add(x)
                rels.append(('merged', tuple(sorted(contributors[x])), x))
        else:
            rels.append(('split', s, tuple(out)))
    for parents, t, role in ((g[1], g[2], g[3]) for g in gen1):
        for x in flow2.get(t, []):
            rels.append(('generated', parents, x, role))
    for g in gen2:
        assert all(p[0] != 'entity' or p[1] not in {r[2] for r in gen1} for p in g[1])
        rels.append(g)
    return sorted(rels, key=relation_sort_key)


def _unique(bodies):
    ids = [e.id for b in bodies for e in b.entities]+[b.id for b in bodies]
    assert len(ids) == len(set(ids)), 'id collision'
