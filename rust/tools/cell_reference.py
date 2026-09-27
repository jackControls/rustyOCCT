"""Independent reference for the cell-complex topology model (TOPOLOGY_MODEL.md).

The seamed `brep_reference.Model` stays the source of OCCT rows. `to_cell`
converts it by rule into the seamless neutral model: a seam pair on a
cylinder whose uses lie a period apart and continue their neighbours in UV is
merged, so each run of the loop between seam uses becomes its own loop with
an integer winding number; seam vertices disappear and full circles without
vertices become ring edges. Any other seam is kept, and the validator rejects
it. Old shell k becomes shell k of the solid region (listing front sides);
its opposite sides form a twin shell: the infinite void's for shell 0, a
bounded void region's for each cavity.

`validate(cell)` is the independent implementation of the cell-model
contract; the Rust validator must reproduce its complete sorted issue lists.
"""
from dataclasses import dataclass, field
import copy
import functools

import mpmath as mp

from fractions import Fraction as F

import spline_cell_reference as spline
from spline_cell_reference import BSpline2, BSpline3, BSplineSurface
from brep_reference import (Arc2, Arc3, cos_rn, sin_rn, Cone, Cylinder, Line2, Line3, Plane, Sphere, TAU, Torus,
                            add, apex, periodic_v,
                            apex_v, axes, cross, periodic, u_scale, v_scale,
                            curve_point, curve_valid, deviation_bounds, dot, finite, margin_check,
                            mul, norm, number, pcurve_point, pcurve_valid, sub, surface_point,
                            surface_valid, vec)


@dataclass
class CEdge:
    start: object      # vertex index or None (ring edge)
    end: object
    curve: object
    fins: list = field(default_factory=list)


@dataclass
class Fin:
    edge: int
    forward: bool
    pcurve: object


@dataclass
class Loop:
    fins: list = field(default_factory=list)   # fin ids; empty with vertex set means a vertex loop
    winding: int = 0                           # turns in u on a cylinder
    vertex: object = None
    winding_v: int = 0                         # turns in v (a torus)


@dataclass
class CFace:
    surface: object
    forward: bool
    loops: list
    front: int
    back: int


@dataclass
class Shell:
    region: int
    sides: list          # (face, 'F' | 'B')
    wire_edges: list = field(default_factory=list)
    acorns: list = field(default_factory=list)


@dataclass
class Region:
    kind: str            # 'solid' | 'void'
    shells: list


@dataclass
class Cell:
    name: str
    tolerance: float
    vertices: list = field(default_factory=list)
    edges: list = field(default_factory=list)
    fins: list = field(default_factory=list)
    loops: list = field(default_factory=list)
    faces: list = field(default_factory=list)
    shells: list = field(default_factory=list)
    regions: list = field(default_factory=list)
    # Declared enclosures (M5): ('v', i), ('u', fin) or ('f', i) to a bound,
    # or to None for a missing one; `declare` fills every absent key.
    enclosures: dict = field(default_factory=dict)


# ---------------------------------------------------------------- enclosures

def gap_bounds(c):
    """{key: (low, high)} of every vertex, fin and face gap, as `validate`
    measures them: vertex to curve ends and vertex-loop surfaces, each use's
    deviation (sampled lower bound, harmonic upper bound or inf where it is
    not harmonic), and consecutive fins' UV gaps (angles scaled by the
    radius). Vertex and UV gaps are computed directly, so low equals high."""
    zero = (mp.mpf(0), mp.mpf(0))
    out = {('v', i): zero for i in range(len(c.vertices))}
    out.update({('f', i): zero for i in range(len(c.faces))})

    def raise_to(key, value):
        try:
            v = value()
            low, high = v if isinstance(v, tuple) else (v, v)
        except (IndexError, TypeError, ZeroDivisionError, ValueError, AttributeError):
            low, high = mp.mpf(0), mp.inf
        old = out.get(key, zero)
        out[key] = (max(old[0], low), max(old[1], high))
    for e in c.edges:
        for v, t in ((e.start, 0), (e.end, 1)):
            if v is not None and 0 <= v < len(c.vertices):
                raise_to(('v', v), lambda: norm(sub(curve_point(e.curve, t), vec(c.vertices[v]))))
    for fi, f in enumerate(c.faces):
        for lid in f.loops:
            if not 0 <= lid < len(c.loops):
                continue
            loop = c.loops[lid]
            if loop.vertex is not None:
                if 0 <= loop.vertex < len(c.vertices):
                    raise_to(('v', loop.vertex), lambda: surface_distance(f.surface, vec(c.vertices[loop.vertex])))
                    if pole_loop(c, f) == lid:
                        raise_to(('v', loop.vertex), lambda: norm(sub(vec(c.vertices[loop.vertex]), pole_point(c, f))))
                continue
            for ui, k in enumerate(loop.fins):
                if not 0 <= k < len(c.fins):
                    continue
                u = c.fins[k]
                raise_to(('u', k), lambda: use_bounds(c.edges[u.edge].curve, f.surface, u.pcurve, u.forward,
                                                      c.tolerance))

                def gap():
                    w = loop.fins[(ui+1) % len(loop.fins)]
                    if isinstance(f.surface, BSplineSurface):
                        d = spline_gap(f.surface, u.pcurve, c.fins[w].pcurve)
                        if d is None:
                            raise ValueError('no exact gap')
                        return d
                    a, b = pcurve_point(u.pcurve, 1), pcurve_point(c.fins[w].pcurve, 0)
                    last = ui == len(loop.fins)-1
                    shift = TAU*loop.winding if last and periodic(f.surface) else 0
                    shift_v = TAU*loop.winding_v if last and periodic(f.surface) else 0
                    du, dv = a[0]-b[0]-shift, a[1]-b[1]-shift_v
                    if periodic(f.surface):
                        du *= u_scale(f.surface, a[1])
                        dv *= v_scale(f.surface)
                    return mp.sqrt(du*du+dv*dv)
                raise_to(('f', fi), gap)
    return out


def use_bounds(curve, s, p, forward, tol):
    """A use's deviation bounds: exact for a rational spline use, by Taylor
    bounds for another."""
    if spline_use(curve, s, p):
        bounds = spline_deviation(curve, s, p, forward)
        if bounds is None:
            bounds = taylor_deviation(curve, s, p, forward, tol)
        if bounds is not None:
            return bounds
    return deviation_bounds(curve, s, p, forward)


def declare(c):
    """Declare every absent enclosure as a generous, sound bound: twice the
    highest gap plus 2^-20 of the tolerance, rounded up and capped at the
    tolerance. (A bound the checker cannot verify is what the mutation cases
    set explicitly.)"""
    import math
    for key, (_, high) in gap_bounds(c).items():
        if key in c.enclosures:
            continue
        if not mp.isfinite(high):
            c.enclosures[key] = c.tolerance
            continue
        bound = math.nextafter(float(2*high+mp.mpf(c.tolerance)*mp.mpf(2)**-20), math.inf)
        c.enclosures[key] = min(bound, c.tolerance)
    return c


def pole_loop(c, f):
    """The loop id of a cone's or sphere's pole: its first vertex loop, when
    its edge loops wind once in total (the pole closes the band)."""
    if not isinstance(f.surface, (Cone, Sphere)):
        return None
    if abs(face_turns(c, f)) != 1:
        return None
    return next((l for l in f.loops if 0 <= l < len(c.loops) and c.loops[l].vertex is not None), None)


def face_turns(c, f):
    return sum(c.loops[l].winding for l in f.loops if 0 <= l < len(c.loops) and c.loops[l].vertex is None)


def pole_north(c, f):
    """The material lies left of a loop's traversal: a band winding +u on a
    forward face closes at its larger v."""
    return (face_turns(c, f) > 0) == f.forward


def pole_v(c, f):
    if isinstance(f.surface, Sphere):
        return mp.pi/2 if pole_north(c, f) else -mp.pi/2
    return apex_v(f.surface)


def pole_point(c, f):
    s = f.surface
    if isinstance(s, Sphere):
        o, x, y, n = axes(s.frame)
        return add(o, mul(n, mp.mpf(s.radius)*(1 if pole_north(c, f) else -1)))
    return apex(s)


# ---------------------------------------------------------------- conversion

def closed_curve(c):
    if isinstance(c, BSpline3):
        return spline.closed_period(c)  # a spline ring edge is a full period
    return isinstance(c, Arc3) and abs(c.sweep) == TAU


def spline_use(curve, surface, pcurve):
    """Spline geometry that only the endpoint and continuity checks
    certify before the rest of S4."""
    return isinstance(curve, BSpline3) or isinstance(surface, BSplineSurface) or isinstance(pcurve, BSpline2)


def spline_gap(s, p, q):
    """The 3D distance between a spline surface's exact points at p's end
    and q's start, or None when an end is not exact or off the domain."""
    a, b = exact_end(p, 1), exact_end(q, 0)
    if a is None or b is None:
        return None
    try:
        x, y = spline.surface_jet(s, *a)[0], spline.surface_jet(s, *b)[0]
    except ValueError:
        return None
    return mp.sqrt(mp_of(sum((i-j)**2 for i, j in zip(x, y))))


def exact_end(p, t):
    """A pcurve's end as Fractions, or None when it is not exact (an arc)."""
    if isinstance(p, BSpline2):
        return spline.end_point(p, t)
    if isinstance(p, Line2):
        return tuple(F(x) for x in (p.start if t == 0 else p.end))
    return None


def seam_merge(face, loop, edges):
    """Split one seamed loop into runs; None when no consistent seam pair."""
    if not isinstance(face.surface, Cylinder):
        return None
    counts = {}
    for u in loop:
        counts[u.edge] = counts.get(u.edge, 0)+1
    seams = [e for e, n in counts.items() if n == 2]
    if not seams:
        return None
    n = len(loop)
    for e in seams:
        a, b = [k for k, u in enumerate(loop) if u.edge == e]
        ua, ub = loop[a], loop[b]
        pa, pb = ua.pcurve, ub.pcurve
        if ua.forward == ub.forward or not isinstance(pa, Line2) or not isinstance(pb, Line2):
            return None
        shift = pa.start[0]-pb.end[0]
        if abs(shift) != TAU or pa.end[0]-pb.start[0] != shift:
            return None
        if pa.start[1] != pb.end[1] or pa.end[1] != pb.start[1]:
            return None
        for k in (a, b):
            prev, nxt = loop[(k-1) % n], loop[(k+1) % n]
            if pcurve_end(prev.pcurve) != loop[k].pcurve.start or loop[k].pcurve.end != pcurve_start(nxt.pcurve):
                return None
    keep = [k for k, u in enumerate(loop) if u.edge not in seams]
    if not keep:
        return None
    # Runs of consecutive kept uses, in loop order starting after a seam use.
    first = next(k for k in range(n) if loop[k].edge in seams)
    runs, run = [], []
    for step in range(1, n+1):
        k = (first+step) % n
        if loop[k].edge in seams:
            if run:
                runs.append(run)
            run = []
        else:
            run.append(k)
    if run:
        runs.append(run)
    out = []
    for run in runs:
        du = pcurve_end(loop[run[-1]].pcurve)[0]-pcurve_start(loop[run[0]].pcurve)[0]
        w = du/TAU
        if w != round(w) or round(w) == 0:
            return None
        out.append(([loop[k] for k in run], int(round(w))))
    return seams, out


def pcurve_start(p):
    if isinstance(p, Line2):
        return p.start
    return (p.center[0]+p.radius*cos_rn(p.start), p.center[1]+p.radius*sin_rn(p.start))


def pcurve_end(p):
    if isinstance(p, Line2):
        return p.end
    a = p.start+p.sweep
    return (p.center[0]+p.radius*cos_rn(a), p.center[1]+p.radius*sin_rn(a))


def to_cell(m):
    """The seamless neutral model of a seamed model, by rule."""
    # Loops per face, with seam pairs merged where consistent.
    face_loops, removed = [], set()
    for f in m.faces:
        loops = []
        for loop in f.loops:
            merged = seam_merge(f, loop, m.edges) if loop else None
            if merged is None:
                loops.append((list(loop), 0))
            else:
                seams, runs = merged
                removed |= set(seams)
                loops.extend(runs)
        face_loops.append(loops)
    # Only edges that every use dropped disappear.
    still_used = {u.edge for loops in face_loops for loop, _ in loops for u in loop}
    removed -= still_used
    edge_map, edges = {}, []
    for i, e in enumerate(m.edges):
        if i not in removed:
            edge_map[i] = len(edges)
            edges.append(copy.deepcopy(e))
    # Vertices used only by removed seams and by closed curves disappear.
    other_use = set()
    for i, e in enumerate(m.edges):
        if i in removed:
            continue
        if not (closed_curve(e.curve) and e.start == e.end):
            other_use |= {e.start, e.end}
    seam_vertices = {v for i in removed for v in (m.edges[i].start, m.edges[i].end)}
    drop = {v for v in seam_vertices if v not in other_use}
    vertex_map, vertices = {}, []
    for v, p in enumerate(m.vertices):
        if v not in drop:
            vertex_map[v] = len(vertices)
            vertices.append(p)
    cell = Cell(m.name, m.tolerance, vertices)
    for e in edges:
        ring = closed_curve(e.curve) and e.start in drop
        start = None if ring else vertex_map.get(e.start, e.start)
        end = None if ring else vertex_map.get(e.end, e.end)
        cell.edges.append(CEdge(start, end, e.curve))
    # Shells: old shell k is shell k of the solid region (front sides); the
    # opposite sides of each non-empty old shell form its twin, after them.
    old = m.shells
    n = len(old)
    twins = {}
    for k, s in enumerate(old):
        if s:
            twins[k] = n+len(twins)
    owner = {}
    for k, s in enumerate(old):
        for f in s:
            owner.setdefault(f, k)
    for fi, f in enumerate(m.faces):
        k = owner.get(fi, 0)
        loop_ids = []
        for loop, w in face_loops[fi]:
            fins = []
            for u in loop:
                fins.append(len(cell.fins))
                cell.fins.append(Fin(edge_map.get(u.edge, u.edge), u.forward, u.pcurve))
            loop_ids.append(len(cell.loops))
            cell.loops.append(Loop(fins, w))
        cell.faces.append(CFace(f.surface, f.forward, loop_ids, k, twins.get(k, k)))
    for fid, fin in enumerate(cell.fins):
        if 0 <= fin.edge < len(cell.edges):
            cell.edges[fin.edge].fins.append(fid)
    kind = getattr(m, 'kind', 'solid')
    if kind == 'sheet':
        # S6: both sides of every face in one shell of the infinite void.
        for f in cell.faces:
            f.front = f.back = 0
        cell.shells = [Shell(0, [(fi, side) for fi in range(len(cell.faces)) for side in 'FB'])]
        cell.regions = [Region('void', [0])]
        return cell
    if kind in ('wire', 'acorn'):
        cell.shells = [Shell(0, [], [edge_map.get(e, e) for e in m.wire] if kind == 'wire' else [],
                             [0] if kind == 'acorn' else [])]
        cell.regions = [Region('void', [0])]
        return cell
    # Region 0 is the infinite void (the outer shell's twin), region 1 the
    # solid (a bounded void for a closed shell without one, S6), then one
    # bounded void region per non-empty cavity.
    cell.regions = [Region('void', [twins[0]] if 0 in twins else []),
                    Region('void' if kind == 'shell' else 'solid', list(range(n)))]
    cell.shells = [Shell(1, [(f, 'F') for f in s]) for s in old]
    for k, t in sorted(twins.items(), key=lambda kv: kv[1]):
        if k == 0:
            region = 0
        else:
            region = len(cell.regions)
            cell.regions.append(Region('void', [t]))
        valid = [f for f in old[k] if 0 <= f < len(m.faces)]
        cell.shells.append(Shell(region, [(f, 'B') for f in valid]))
    return cell


# ---------------------------------------------------------------- protocol

def encode(c):
    """Line protocol read by the Rust fixture test."""
    out = [f'case {c.name}', f'tolerance {number(c.tolerance)}']
    frame = lambda f: ' '.join(map(number, (*f.origin, *f.normal, *f.x)))
    ref = lambda v: '-' if v is None else str(v)
    enc = lambda key: ' enc '+('-' if c.enclosures.get(key) is None else number(c.enclosures[key]))
    for i, v in enumerate(c.vertices):
        out.append('v '+' '.join(map(number, v))+enc(('v', i)))
    for e in c.edges:
        cv = e.curve
        head = f'e {ref(e.start)} {ref(e.end)}'
        if isinstance(cv, Line3):
            body = 'line '+' '.join(map(number, (*cv.start, *cv.end)))
        elif isinstance(cv, BSpline3):
            body = spline.encode_curve(cv, number)
        else:
            body = f'arc {frame(cv.frame)} '+' '.join(map(number, (cv.radius, cv.start, cv.sweep)))
        out.append(f'{head} {body} fins'+''.join(f' {k}' for k in e.fins))
    # Loops (with their fins) in arena order; faces name their loops.
    for loop in c.loops:
        if loop.vertex is not None:
            out.append(f'lv {loop.vertex}')
            continue
        out.append(f'l {loop.winding}' + (f' {loop.winding_v}' if loop.winding_v else ''))
        for k in loop.fins:
            u = c.fins[k]
            p = u.pcurve
            o = 'F' if u.forward else 'R'
            if isinstance(p, Line2):
                out.append(f'u {u.edge} {o} line '+' '.join(map(number, (*p.start, *p.end)))+enc(('u', k)))
            elif isinstance(p, BSpline2):
                out.append(f'u {u.edge} {o} {spline.encode_curve(p, number)}'+enc(('u', k)))
            else:
                out.append(f'u {u.edge} {o} arc '
                           + ' '.join(map(number, (*p.center, p.radius, p.start, p.sweep)))+enc(('u', k)))
    for fi, f in enumerate(c.faces):
        s = f.surface
        o = 'F' if f.forward else 'R'
        loops = ' loops'+''.join(f' {l}' for l in f.loops)
        if isinstance(s, Plane):
            out.append(f'f plane {frame(s.frame)} {o} {f.front} {f.back}{loops}'+enc(('f', fi)))
        elif isinstance(s, BSplineSurface):
            out.append(f'f {spline.encode_surface(s, number)} {o} {f.front} {f.back}{loops}'+enc(('f', fi)))
        elif isinstance(s, Cone):
            out.append(f'f cone {frame(s.frame)} {number(s.radius)} {number(s.half_angle)} {o} {f.front} {f.back}'
                       f'{loops}'+enc(('f', fi)))
        elif isinstance(s, Sphere):
            out.append(f'f sphere {frame(s.frame)} {number(s.radius)} {o} {f.front} {f.back}{loops}'
                       + enc(('f', fi)))
        elif isinstance(s, Torus):
            out.append(f'f torus {frame(s.frame)} {number(s.major)} {number(s.minor)} {o} {f.front} {f.back}'
                       f'{loops}'+enc(('f', fi)))
        else:
            out.append(f'f cylinder {frame(s.frame)} {number(s.radius)} {o} {f.front} {f.back}{loops}'
                       + enc(('f', fi)))
    for s in c.shells:
        sides = ' '.join(f'{f}:{side}' for f, side in s.sides)
        out.append(f's {s.region} sides {sides} wire'+''.join(f' {e}' for e in s.wire_edges)
                   + ' acorn'+''.join(f' {v}' for v in s.acorns))
    for r in c.regions:
        out.append(f'r {r.kind}'+''.join(f' {s}' for s in r.shells))
    out.append('end')
    return '\n'.join(out)


# ---------------------------------------------------------------- validator

def body_class(c):
    """D9 (S6): solid, sheet, wire, acorn or general, from the regions and
    what the shells list."""
    wires = any(s.wire_edges for s in c.shells)
    acorns = any(s.acorns for s in c.shells)
    solid = [r for r in c.regions if r.kind == 'solid']
    kind = lambda shell: c.regions[c.shells[shell].region].kind
    if c.faces and not wires and not acorns:
        sides = [(kind(f.front), kind(f.back)) for f in c.faces]
        if solid and all(sorted(p) == ['solid', 'void'] for p in sides):
            return 'solid'
        if not solid:
            return 'sheet'
    if not c.faces and not solid:
        if wires and not acorns:
            return 'wire'
        if acorns and not wires and sum(len(s.acorns) for s in c.shells) == 1:
            return 'acorn'
    return 'general'


def issue(kind, entity):
    return (kind, entity)


def fin_vertices(c, k):
    fin = c.fins[k]
    e = c.edges[fin.edge]
    return (e.start, e.end) if fin.forward else (e.end, e.start)


def opposite(side):
    return 'B' if side == 'F' else 'F'


def validate(c):
    """Complete sorted issue list for the cell-model contract."""
    issues = []
    nv, ne, nfin, nl, nf, ns, nr = (len(c.vertices), len(c.edges), len(c.fins), len(c.loops),
                                    len(c.faces), len(c.shells), len(c.regions))
    rng = lambda i, n: i is not None and 0 <= i < n
    # Fin and loop locations (by first occurrence) name the entities.
    where, loop_face = {}, {}
    for fi, f in enumerate(c.faces):
        for li, lid in enumerate(f.loops):
            if rng(lid, nl):
                loop_face.setdefault(lid, (fi, li))
                for ui, k in enumerate(c.loops[lid].fins):
                    if rng(k, nfin):
                        where.setdefault(k, (fi, li, ui))
    fin_name = lambda k: f'use {where[k][0]}.{where[k][1]}.{where[k][2]}' if k in where else f'fin {k}'
    loop_name = lambda lid: f'loop {loop_face[lid][0]}.{loop_face[lid][1]}' if lid in loop_face else f'loop slot {lid}'

    # ------------------------------------------------ references
    for i, e in enumerate(c.edges):
        if any(v is not None and not rng(v, nv) for v in (e.start, e.end)) or any(not rng(k, nfin) for k in e.fins):
            issues.append(issue('reference', f'edge {i}'))
    for k, fin in enumerate(c.fins):
        if not rng(fin.edge, ne):
            issues.append(issue('reference', fin_name(k)))
    for lid, loop in enumerate(c.loops):
        if any(not rng(k, nfin) for k in loop.fins) or (loop.vertex is not None and not rng(loop.vertex, nv)):
            issues.append(issue('reference', loop_name(lid)))
    for fi, f in enumerate(c.faces):
        if any(not rng(l, nl) for l in f.loops) or not rng(f.front, ns) or not rng(f.back, ns):
            issues.append(issue('reference', f'face {fi}'))
    for si, s in enumerate(c.shells):
        if (not rng(s.region, nr) or any(not rng(f, nf) for f, _ in s.sides)
                or any(not rng(e, ne) for e in s.wire_edges) or any(not rng(v, nv) for v in s.acorns)):
            issues.append(issue('reference', f'shell {si}'))
    for ri, r in enumerate(c.regions):
        if any(not rng(s, ns) for s in r.shells):
            issues.append(issue('reference', f'region {ri}'))
    if issues:
        return sorted(set(issues))
    tol = mp.mpf(c.tolerance)

    # ------------------------------------------------ structure (exact)
    fin_count = {}
    for loop in c.loops:
        for k in loop.fins:
            fin_count[k] = fin_count.get(k, 0)+1
    for k in range(nfin):
        if fin_count.get(k, 0) == 0:
            issues.append(issue('fin_without_loop', fin_name(k)))
        elif fin_count[k] > 1:
            issues.append(issue('fin_reused', fin_name(k)))
    loop_count = {}
    for f in c.faces:
        for lid in f.loops:
            loop_count[lid] = loop_count.get(lid, 0)+1
    for lid in range(nl):
        if loop_count.get(lid, 0) == 0:
            issues.append(issue('loop_without_face', loop_name(lid)))
        elif loop_count[lid] > 1:
            issues.append(issue('loop_reused', loop_name(lid)))
    users = {}
    for k, fin in enumerate(c.fins):
        users.setdefault(fin.edge, []).append(k)
    for i, e in enumerate(c.edges):
        if sorted(e.fins) != sorted(users.get(i, [])):
            issues.append(issue('edge_fins_mismatch', f'edge {i}'))
    used_vertices = {v for e in c.edges for v in (e.start, e.end) if v is not None}
    used_vertices |= {l.vertex for l in c.loops if l.vertex is not None}
    used_vertices |= {v for s in c.shells for v in s.acorns}
    for v in range(nv):
        if v not in used_vertices:
            issues.append(issue('unused_vertex', f'vertex {v}'))
    wire = {e for s in c.shells for e in s.wire_edges}
    for i in range(ne):
        if i not in users and i not in wire:
            issues.append(issue('unused_edge', f'edge {i}'))
        e = c.edges[i]
        if (e.start is None) != (e.end is None):
            issues.append(issue('ring_edge_with_vertex', f'edge {i}'))
        elif e.start is None and not closed_curve(e.curve):
            issues.append(issue('ring_edge_open', f'edge {i}'))
    # Face sides against shells.
    listings = {}
    for si, s in enumerate(c.shells):
        if not s.sides and not s.wire_edges and not s.acorns:
            issues.append(issue('empty_shell', f'shell {si}'))
        for f, side in s.sides:
            listings.setdefault((f, side), []).append(si)
    side_bad = set()
    for fi, f in enumerate(c.faces):
        fr, bk = listings.get((fi, 'F'), []), listings.get((fi, 'B'), [])
        if not fr and not bk:
            issues.append(issue('face_without_shell', f'face {fi}'))
        elif len(fr) > 1 or len(bk) > 1:
            same_shell = any(fr.count(s) > 1 for s in fr) or any(bk.count(s) > 1 for s in bk)
            issues.append(issue('face_reused' if same_shell else 'side_in_two_shells', f'face {fi}'))
        elif not fr or not bk:
            issues.append(issue('side_without_shell', f'face {fi}'))
        elif fr[0] != f.front or bk[0] != f.back:
            issues.append(issue('side_region_mismatch', f'face {fi}'))
        else:
            continue
        side_bad.add(fi)
    for fi, f in enumerate(c.faces):
        # A sphere or torus without loops is the whole closed surface.
        if not f.loops and not isinstance(f.surface, (Sphere, Torus)):
            issues.append(issue('empty_face', f'face {fi}'))
        for li, lid in enumerate(f.loops):
            loop = c.loops[lid]
            name = f'loop {fi}.{li}'
            if loop.vertex is not None:
                continue
            if not loop.fins:
                issues.append(issue('empty_loop', name))
                continue
            if (not periodic(f.surface) and loop.winding != 0) or (not periodic_v(f.surface) and loop.winding_v != 0):
                issues.append(issue('winding_mismatch', name))
            ends = [fin_vertices(c, k) for k in loop.fins]
            if len(ends) == 1 and ends[0] == (None, None):
                continue
            if any(None in e for e in ends) or any(ends[k][1] != ends[(k+1) % len(ends)][0] for k in range(len(ends))):
                issues.append(issue('open_loop', name))
        if periodic(f.surface):
            # Windings balance, except on a cone where one pole (a vertex
            # loop at the apex) closes a band that winds once.
            wound = [c.loops[lid].winding for lid in f.loops if c.loops[lid].vertex is None]
            total = sum(wound)
            total_v = sum(c.loops[lid].winding_v for lid in f.loops if c.loops[lid].vertex is None)
            if (total != 0 and not (isinstance(f.surface, (Cone, Sphere)) and abs(total) == 1
                                    and pole_loop(c, f) is not None)) or total_v != 0:
                issues.append(issue('winding_mismatch', f'loop {fi}.0'))
    for ri, r in enumerate(c.regions):
        if len(set(r.shells)) != len(r.shells):
            issues.append(issue('double_bounding', f'region {ri}'))
        if ri > 0 and not r.shells:
            issues.append(issue('region_without_shell', f'region {ri}'))
    if not c.regions or c.regions[0].kind != 'void':
        issues.append(issue('no_infinite_region', 'region 0'))
    for si, s in enumerate(c.shells):
        if si not in c.regions[s.region].shells:
            issues.append(issue('region_shell_mismatch', f'shell {si}'))
    # Seams are forbidden: an edge with two fins in one face.
    fin_face = {k: where[k][0] for k in where}
    for i in range(ne):
        faces = [fin_face.get(k) for k in users.get(i, [])]
        if len(faces) != len(set(faces)):
            issues.append(issue('seam_edge', f'edge {i}'))

    # Edge accounting: shells must alternate around each edge's fins.
    bad_shells = set()
    for i, e in enumerate(c.edges):
        fins = [k for k in e.fins if k in fin_face and fin_face[k] not in side_bad]
        if not fins or sorted(e.fins) != sorted(users.get(i, [])):
            continue
        faces = [c.faces[fin_face[k]] for k in fins]
        ahead = [f.front if c.fins[k].forward else f.back for k, f in zip(fins, faces)]
        behind = [f.back if c.fins[k].forward else f.front for k, f in zip(fins, faces)]
        n = len(fins)
        if all(ahead[j] == behind[(j+1) % n] for j in range(n)):
            continue
        pairs = {frozenset((f.front, f.back)) for f in faces}
        if len(pairs) > 1:
            kind = 'edge_across_shells'
        elif n == 1:
            kind = 'free_edge'
        elif n == 2:
            kind = 'same_sense_uses' if c.fins[fins[0]].forward == c.fins[fins[1]].forward else 'radial_order_inconsistent'
        else:
            kind = 'non_manifold_edge'
        issues.append(issue(kind, f'edge {i}'))
        for f in faces:
            bad_shells |= {f.front, f.back}
    structural_faces = set()
    for kind, ent in issues:
        if kind in ('open_loop', 'empty_loop', 'empty_face', 'winding_mismatch'):
            structural_faces.add(int(ent.split()[1].split('.')[0]))
    for i in range(ne):
        if any(k == 'seam_edge' and ent == f'edge {i}' for k, ent in issues):
            structural_faces |= {fin_face[k] for k in users.get(i, []) if k in fin_face}
    for fi in structural_faces:
        bad_shells |= {c.faces[fi].front, c.faces[fi].back}

    def twin_of(si):
        mine = sorted((f, opposite(side)) for f, side in c.shells[si].sides)
        for sj in range(si):
            if mine and sorted(c.shells[sj].sides) == mine:
                return sj
        return None

    twins = {si for si in range(ns) if twin_of(si) is not None}
    shell_faces = lambda si: [f for f, _ in c.shells[si].sides]
    for si, s in enumerate(c.shells):
        if si in twins:
            continue
        faces = shell_faces(si)
        if si in bad_shells or not faces or any(f in side_bad for f in faces):
            bad_shells.add(si)
            continue
        members = set(faces)
        adjacency = {f: set() for f in members}
        for i in range(ne):
            fs = {fin_face[k] for k in users.get(i, []) if k in fin_face and fin_face[k] in members}
            for a in fs:
                adjacency[a] |= fs-{a}
        seen, stack = set(), [faces[0]]
        while stack:
            f = stack.pop()
            if f not in seen:
                seen.add(f)
                stack.extend(adjacency[f])
        if seen != members:
            issues.append(issue('disconnected_shell', f'shell {si}'))
            bad_shells.add(si)
            continue
        links = {}
        for f in faces:
            for lid in c.faces[f].loops:
                loop = c.loops[lid]
                for j, k in enumerate(loop.fins):
                    w = loop.fins[(j+1) % len(loop.fins)]
                    v = fin_vertices(c, k)[1]
                    if v is None:
                        continue
                    g = links.setdefault(v, {})
                    g.setdefault(c.fins[k].edge, set()).add(c.fins[w].edge)
                    g.setdefault(c.fins[w].edge, set()).add(c.fins[k].edge)
        pinched = False
        for v in sorted(links):
            g = links[v]
            start = next(iter(g))
            seen, stack = set(), [start]
            while stack:
                x = stack.pop()
                if x not in seen:
                    seen.add(x)
                    stack.extend(g[x])
            if seen != set(g):
                issues.append(issue('non_manifold_vertex', f'vertex {v}'))
                pinched = True
        if pinched:
            bad_shells.add(si)
            continue
        vs, es, loops = set(), set(), 0
        for f in members:
            for lid in c.faces[f].loops:
                loop = c.loops[lid]
                loops += 1
                if loop.vertex is not None:
                    vs.add(loop.vertex)
                for k in loop.fins:
                    edge = c.edges[c.fins[k].edge]
                    if edge.start is not None and edge.end is not None:
                        vs |= {edge.start, edge.end}
                        es.add(c.fins[k].edge)
        chi = len(vs)-len(es)+2*len(members)-loops
        # A sheet's two sides make no closed surface of its faces (S6).
        if not two_sided(c, si) and (chi % 2 or chi > 2):
            issues.append(issue('euler', f'shell {si}'))
            bad_shells.add(si)

    # Wire edges and acorn vertices (S6): a wire edge has no fins, an acorn
    # vertex bounds nothing else, and a shell's wire edges and acorns are
    # connected through their vertices.
    edge_users = {v for e in c.edges for v in (e.start, e.end) if v is not None}
    loop_vertices = {l.vertex for l in c.loops if l.vertex is not None}
    for si, sh in enumerate(c.shells):
        for e in sh.wire_edges:
            if c.edges[e].fins or e in users:
                issues.append(issue('wire_edge_with_fins', f'edge {e}'))
        for v in sh.acorns:
            if v in edge_users or v in loop_vertices:
                issues.append(issue('acorn_vertex_used', f'vertex {v}'))
        # A shell of wire edges and acorns alone is connected through their
        # vertices (a mixed general body is not constrained here).
        if (sh.wire_edges or sh.acorns) and not sh.sides:
            parts = [('e', e) for e in sh.wire_edges]+[('v', v) for v in sh.acorns]
            ends = lambda p: ({c.edges[p[1]].start, c.edges[p[1]].end} - {None}) if p[0] == 'e' else {p[1]}
            seen, stack = {parts[0]}, [parts[0]]
            while stack:
                p = stack.pop()
                for q in parts:
                    if q not in seen and ends(p) & ends(q):
                        seen.add(q)
                        stack.append(q)
            if len(seen) != len(parts):
                issues.append(issue('disconnected_shell', f'shell {si}'))

    # ------------------------------------------------ certified-in-production geometry
    vertex_ok = []
    for v, p in enumerate(c.vertices):
        vertex_ok.append(finite(p))
        if not vertex_ok[-1]:
            issues.append(issue('degenerate_vertex', f'vertex {v}'))
    curve_ok = []
    for i, e in enumerate(c.edges):
        ok = curve_valid(e.curve, tol)
        curve_ok.append(ok)
        if not ok:
            issues.append(issue('degenerate_curve', f'edge {i}'))
    surface_ok = []
    for fi, f in enumerate(c.faces):
        ok = surface_valid(f.surface, tol)
        surface_ok.append(ok)
        if not ok:
            issues.append(issue('degenerate_surface', f'face {fi}'))
    geometry_bad = {fi for fi, ok in enumerate(surface_ok) if not ok}
    # C1 in each cell's own parameterisation (R4), exactly.
    for i, e in enumerate(c.edges):
        if isinstance(e.curve, BSpline3) and not spline.curve_c1(e.curve):
            issues.append(issue('edge_not_c1', f'edge {i}'))
    for fi, f in enumerate(c.faces):
        if isinstance(f.surface, BSplineSurface) and not spline.surface_c1(f.surface):
            issues.append(issue('face_not_c1', f'face {fi}'))
            geometry_bad.add(fi)
        for li, lid in enumerate(f.loops):
            for ui, k in enumerate(c.loops[lid].fins):
                p = c.fins[k].pcurve
                if isinstance(p, BSpline2) and not spline.curve_c1(p):
                    issues.append(issue('pcurve_not_c1', f'use {fi}.{li}.{ui}'))
                    geometry_bad.add(fi)
    # Enclosures (M5): a usable bound lies in [0, tol]; each geometric check
    # below decides against it first, then against the tolerance.
    def usable(key, entity):
        if c.enclosures.get(key) is None:
            issues.append(issue('enclosure_missing', entity))
            return None
        b = c.enclosures[key]
        if not (0 <= b <= tol):
            issues.append(issue('enclosure_exceeds_resolution', entity))
            return None
        return mp.mpf(b)
    vertex_bound = [usable(('v', v), f'vertex {v}') for v in range(nv)]
    face_bound = [usable(('f', fi), f'face {fi}') for fi in range(nf)]
    fin_bound = {k: usable(('u', k), fin_name(k)) for k in sorted(where, key=lambda k: where[k])}

    def judge(low, high, bound, entity, what):
        """'within' or 'beyond' the tolerance, reporting the bound's verdict."""
        if bound is not None:
            margin_check(low, bound, what)
            margin_check(high, bound, what)
            if high <= bound:
                return 'within'
        margin_check(low, tol, what)
        margin_check(high, tol, what)
        if low > tol:
            return 'beyond'
        if high > tol:
            raise ArithmeticError(f'{c.name}: {what} is too close to tolerance for the oracle')
        if bound is not None:
            if low > bound:
                issues.append(issue('enclosure_unsound', entity))
            else:
                raise ArithmeticError(f'{c.name}: {what} is too close to its enclosure for the oracle')
        return 'within'

    for i, e in enumerate(c.edges):
        if not curve_ok[i] or e.start is None or e.end is None:
            continue
        for end, (v, t) in enumerate(((e.start, 0), (e.end, 1))):
            if not vertex_ok[v]:
                continue
            d = norm(sub(curve_point(e.curve, t), vec(c.vertices[v])))
            if judge(d, d, vertex_bound[v], f'vertex {v}', f'vertex on edge {i}') == 'beyond':
                issues.append(issue('vertex_off_curve', f'edge {i} {"start" if end == 0 else "end"}'))
    for fi, f in enumerate(c.faces):
        for li, lid in enumerate(f.loops):
            loop = c.loops[lid]
            if loop.vertex is not None:
                if surface_ok[fi] and vertex_ok[loop.vertex] and isinstance(f.surface, BSplineSurface):
                    issues.append(issue('uncertified_vertex_loop', f'loop {fi}.{li}'))
                    geometry_bad.add(fi)
                elif surface_ok[fi] and vertex_ok[loop.vertex]:
                    d = surface_distance(f.surface, vec(c.vertices[loop.vertex]))
                    v = loop.vertex
                    if judge(d, d, vertex_bound[v], f'vertex {v}', f'{c.name}: vertex loop') == 'beyond':
                        issues.append(issue('vertex_loop_off_surface', f'loop {fi}.{li}'))
                        geometry_bad.add(fi)
                    elif pole_loop(c, f) == lid:
                        d = norm(sub(vec(c.vertices[v]), pole_point(c, f)))
                        if judge(d, d, vertex_bound[v], f'vertex {v}', f'{c.name}: pole') == 'beyond':
                            issues.append(issue('pole_off_apex', f'loop {fi}.{li}'))
                            geometry_bad.add(fi)
                continue
            for ui, k in enumerate(loop.fins):
                u = c.fins[k]
                ent = f'use {fi}.{li}.{ui}'
                if not pcurve_valid(u.pcurve):
                    issues.append(issue('degenerate_pcurve', ent))
                    geometry_bad.add(fi)
                    continue
                if not surface_ok[fi] or not curve_ok[u.edge]:
                    continue
                if spline_use(c.edges[u.edge].curve, f.surface, u.pcurve):
                    bounds = spline_deviation(c.edges[u.edge].curve, f.surface, u.pcurve, u.forward)
                    if bounds is None:
                        bounds = taylor_deviation(c.edges[u.edge].curve, f.surface, u.pcurve, u.forward,
                                                  c.tolerance)
                    if bounds is None:
                        issues.append(issue('uncertified_pcurve_off_edge', ent))
                        geometry_bad.add(fi)
                        continue
                    low, high = bounds
                else:
                    low, high = deviation_bounds(c.edges[u.edge].curve, f.surface, u.pcurve, u.forward)
                if judge(low, high, fin_bound.get(k), ent, f'{c.name}: {ent}') == 'beyond':
                    issues.append(issue('pcurve_off_edge', ent))
                    geometry_bad.add(fi)

    # UV continuity; the last fin closes on the first shifted by the winding.
    for fi, f in enumerate(c.faces):
        if not surface_ok[fi]:
            continue
        for li, lid in enumerate(f.loops):
            loop = c.loops[lid]
            if loop.vertex is not None or any(not pcurve_valid(c.fins[k].pcurve) for k in loop.fins):
                continue
            for ui, k in enumerate(loop.fins):
                last = ui == len(loop.fins)-1
                w = loop.fins[(ui+1) % len(loop.fins)]
                if isinstance(f.surface, BSplineSurface):
                    # Measured in 3D between the exact surface points.
                    d = spline_gap(f.surface, c.fins[k].pcurve, c.fins[w].pcurve)
                    if d is None:
                        issues.append(issue('uncertified_uv_gap', f'use {fi}.{li}.{ui}'))
                        geometry_bad.add(fi)
                    elif judge(d, d, face_bound[fi], f'face {fi}', f'{c.name}: uv gap') == 'beyond':
                        issues.append(issue('uv_gap', f'use {fi}.{li}.{ui}'))
                        geometry_bad.add(fi)
                    continue
                a, b = pcurve_point(c.fins[k].pcurve, 1), pcurve_point(c.fins[w].pcurve, 0)
                shift = TAU*loop.winding if last and periodic(f.surface) else 0
                shift_v = TAU*loop.winding_v if last and periodic(f.surface) else 0
                du, dv = a[0]-b[0]-shift, a[1]-b[1]-shift_v
                if periodic(f.surface):
                    du *= u_scale(f.surface, a[1])
                    dv *= v_scale(f.surface)
                d = mp.sqrt(du*du+dv*dv)
                if judge(d, d, face_bound[fi], f'face {fi}', f'{c.name}: uv gap') == 'beyond':
                    issues.append(issue('uv_gap', f'use {fi}.{li}.{ui}'))
                    geometry_bad.add(fi)

    # Loop winding and imbrication, for structurally and geometrically sound faces.
    for fi, f in enumerate(c.faces):
        if fi in geometry_bad or fi in structural_faces or not f.loops:
            continue
        loops = [c.loops[lid] for lid in f.loops]
        sense = 1 if f.forward else -1
        wound = periodic(f.surface) and any(l.winding for l in loops)
        wound_v = periodic_v(f.surface) and any(l.winding_v for l in loops)
        if wound and wound_v:
            issues.append(issue('uncertified_loop_winding', f'loop {fi}.0'))
            continue
        if wound_v:
            total = sum(periodic_area_v(c, l) for l in loops if l.vertex is None)
            if total*sense <= 0:
                issues.append(issue('loop_winding', f'loop {fi}.0'))
            for li, l in enumerate(loops):
                if l.winding_v == 0 and l.vertex is None:
                    if periodic_area_v(c, l)*sense >= 0:
                        issues.append(issue('loop_winding', f'loop {fi}.{li}'))
                    else:
                        issues.append(issue('uncertified_containment', f'loop {fi}.{li}'))
            continue
        if wound:
            total = sum(periodic_area(c, l) for l in loops)
            if pole_loop(c, f) is not None:
                # The pole is the line v = v_apex traversed against the band.
                total += 2*mp.pi*sum(l.winding for l in loops if l.vertex is None)*pole_v(c, f)
            if total*sense <= 0:
                issues.append(issue('loop_winding', f'loop {fi}.0'))
            for li, l in enumerate(loops):
                if l.winding == 0 and l.vertex is None:
                    if periodic_area(c, l)*sense >= 0:
                        issues.append(issue('loop_winding', f'loop {fi}.{li}'))
                    else:
                        issues.append(issue('uncertified_containment', f'loop {fi}.{li}'))
            continue
        # The first edge loop is the outer loop.
        edge_loops = [(li, l) for li, l in enumerate(loops) if l.vertex is None]
        for pos, (li, l) in enumerate(edge_loops):
            a = loop_area(c, l)
            want = (1 if pos == 0 else -1)*sense
            if abs(a) < mp.mpf('1e-25'):
                raise ArithmeticError(f'{c.name}: loop area too small for the oracle')
            if a*want <= 0:
                issues.append(issue('loop_winding', f'loop {fi}.{li}'))
        for li, l in edge_loops[1:]:
            point = pcurve_point(c.fins[l.fins[0]].pcurve, 0)
            if winding(c, edge_loops[0][1], point) == 0:
                issues.append(issue('inner_loop_outside', f'loop {fi}.{li}'))

    # Region orientation and cavity nesting, for fully sound shells.
    bad_faces = (geometry_bad | structural_faces | side_bad
                 | {int(ent.split()[1].split('.')[0]) for k, ent in issues if k == 'loop_winding'})
    sound = {si for si in range(ns) if si not in bad_shells and si not in twins
             and shell_faces(si) and not (set(shell_faces(si)) & bad_faces)
             and set(shell_faces(si)) - two_sided(c, si)}
    for ri, r in enumerate(c.regions):
        oriented = []
        for pos, si in enumerate(r.shells):
            if si not in sound:
                continue
            flux = shell_flux(c, si)
            if flux is None:
                issues.append(issue('uncertified_shell_orientation', f'shell {si}'))
                continue
            if abs(flux) < mp.mpf('1e-25'):
                raise ArithmeticError(f'{c.name}: shell volume too small for the oracle')
            outer = ri > 0 and pos == 0
            if (flux > 0) != outer:
                issues.append(issue('shell_orientation', f'shell {si}'))
            else:
                oriented.append(si)
        if ri == 0 or not r.shells or r.shells[0] not in oriented:
            continue
        outer = r.shells[0]
        cavities = [si for si in oriented if si != outer]
        for si in cavities:
            point = shell_point(c, si, vertex_ok)
            if point is None:
                continue
            # Rays against cones are not decided yet (None): uncertified.
            outside = inside(c, outer, point)
            if outside is None:
                issues.append(issue('uncertified_containment', f'shell {si}'))
                continue
            if not outside:
                issues.append(issue('cavity_outside', f'shell {si}'))
                continue
            nested = False
            for sj in cavities:
                if sj == si:
                    continue
                found = inside(c, sj, point)
                if found is None:
                    issues.append(issue('uncertified_containment', f'shell {si}'))
                    nested = False
                    break
                nested = nested or found
            if nested:
                issues.append(issue('nested_cavity', f'shell {si}'))
    return sorted(set(issues))


def mp_of(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, F) else mp.mpf(x)


def spline_cuts(curve, s, p, forward):
    """The common pieces of a use in the pcurve's fraction: the ends, both
    curves' knots, and a line pcurve's crossings of a spline surface's knot
    lines (the kernel's exact pieces)."""
    cuts = {F(0), F(1)}
    if isinstance(p, BSpline2):
        cuts |= set(spline.knot_fractions(p))
    if isinstance(curve, BSpline3):
        cuts |= {k if forward else 1-k for k in spline.knot_fractions(curve)}
    if isinstance(s, BSplineSurface) and isinstance(p, Line2):
        a, b = [F(x) for x in p.start], [F(x) for x in p.end]
        for box in spline.patch_boxes(s):
            for axis in range(2):
                d = b[axis]-a[axis]
                if d:
                    for k in box[axis]:
                        t = (k-a[axis])/d
                        if 0 < t < 1:
                            cuts.add(t)
    return sorted(cuts)


def pcurve_controls(p, t0, t1):
    """The pcurve's Bézier control points (u, v) on [t0, t1]."""
    if isinstance(p, Line2):
        a, b = [F(x) for x in p.start], [F(x) for x in p.end]
        return [tuple(x+t*(y-x) for x, y in zip(a, b)) for t in (t0, t1)]
    return [(q[0]/q[-1], q[1]/q[-1]) for q in spline.blossom_piece(p, t0, t1)]


def box_of(s, points):
    for box in spline.patch_boxes(s):
        (u0, u1), (v0, v1) = box
        if all(u0 <= u <= u1 and v0 <= v <= v1 for u, v in points):
            return box
    return None


def degree_of(curve):
    return curve.basis.degree if isinstance(curve, (BSpline2, BSpline3)) else 1


def exact_point(curve, t):
    """A line's or spline's point at fraction t, exactly."""
    if isinstance(curve, (Line2, Line3)):
        return tuple(F(x)+fraction_of(t)*(F(y)-F(x)) for x, y in zip(curve.start, curve.end))
    return spline.curve_jet(curve, t)[0]


def fraction_of(t):
    return spline.fraction(t)


def spline_deviation(curve, s, p, forward):
    """(low, high) of a rational spline use, or None where the kernel does
    not decide it: exact samples on every common piece, enough to fix the
    numerator's degree, so zero at all of them is zero everywhere (high 0);
    otherwise the largest sample is certain and high is unbounded."""
    if not (isinstance(curve, (Line3, BSpline3)) and isinstance(p, (Line2, BSpline2))
            and isinstance(s, (Plane, BSplineSurface))):
        return None
    cuts = spline_cuts(curve, s, p, forward)
    if isinstance(s, Plane):
        o, x, y, _ = axes(s.frame)
        surface = lambda u, v: add(o, add(mul(x, u), mul(y, v)))
        extra = 1
    else:
        surface = lambda u, v: [mp_of(z) for z in spline.surface_jet(s, u, v)[0]]
        extra = s.u.degree+s.v.degree
    low = mp.mpf(0)
    for t0, t1 in zip(cuts, cuts[1:]):
        if isinstance(s, BSplineSurface):
            if box_of(s, pcurve_controls(p, t0, t1)) is None:
                return None
            if extra*degree_of(p)+degree_of(curve)+1 > 96:
                return None
        count = extra*degree_of(p)+degree_of(curve)+2
        for k in range(count+1):
            t = t0+(t1-t0)*F(k, count)
            q = exact_point(p, t)
            e = exact_point(curve, t if forward else 1-t)
            d = sub([mp_of(z) for z in e], surface(q[0], q[1]))
            low = max(low, norm(d))
    return (low, mp.mpf(0)) if low < mp.mpf('1e-25') else (low, mp.inf)


def _point_jet(curve, t):
    """(point, derivative) of a line, arc or spline at fraction t (mp)."""
    t = mp.mpf(t)
    if isinstance(curve, (Line2, Line3)):
        a, b = vec(curve.start), vec(curve.end)
        return add(a, mul(sub(b, a), t)), sub(b, a)
    if isinstance(curve, (BSpline2, BSpline3)):
        q, d = spline.curve_jet_mp(curve, t)
        return list(q), list(d)
    theta = mp.mpf(curve.start)+mp.mpf(curve.sweep)*t
    sw, r = mp.mpf(curve.sweep), mp.mpf(curve.radius)
    if isinstance(curve, Arc3):
        o, x, y, _ = axes(curve.frame)
    else:
        o, x, y = vec(curve.center), [mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]
    return (add(o, add(mul(x, r*mp.cos(theta)), mul(y, r*mp.sin(theta)))),
            add(mul(x, -r*sw*mp.sin(theta)), mul(y, r*sw*mp.cos(theta))))


def _interval_jet(curve, t):
    """(value, first, second derivative) of a curve over the interval t."""
    iv = mp.iv
    if isinstance(curve, (Line2, Line3)):
        a, b = [iv.mpf(x) for x in curve.start], [iv.mpf(x) for x in curve.end]
        return [x+(y-x)*t for x, y in zip(a, b)], [y-x for x, y in zip(a, b)], [iv.mpf(0)]*len(a)
    if isinstance(curve, (BSpline2, BSpline3)):
        return spline.spline_jet2_iv(curve, t)
    theta = iv.mpf(curve.start)+iv.mpf(curve.sweep)*t
    sw, r = iv.mpf(curve.sweep), iv.mpf(curve.radius)
    if isinstance(curve, Arc3):
        o, x, y, _ = axes(curve.frame)
        o, x, y = [[iv.mpf(z) for z in v] for v in (o, x, y)]
    else:
        o, x, y = [iv.mpf(z) for z in curve.center], [iv.mpf(1), iv.mpf(0)], [iv.mpf(0), iv.mpf(1)]
    co, si = iv.cos(theta), iv.sin(theta)
    return ([a+r*(b*co+c*si) for a, b, c in zip(o, x, y)],
            [r*sw*(-b*si+c*co) for b, c in zip(x, y)],
            [-r*sw*sw*(b*co+c*si) for b, c in zip(x, y)])


def _surface_jet(s, u, v, interval):
    """S, S_u, S_v (point, mp) or all six partials over a box (intervals)."""
    if isinstance(s, BSplineSurface):
        if interval:
            return spline.surface_jet2_iv(s, u, v)
        return [list(x) for x in spline.surface_jet_mp(s, u, v)]
    ar = mp.iv if interval else mp
    o, x, y, n = [[ar.mpf(z) for z in w] for w in axes(s.frame)]
    co, si = ar.cos(u), ar.sin(u)
    e = [a*co+b*si for a, b in zip(x, y)]
    e1 = [b*co-a*si for a, b in zip(x, y)]
    comb = lambda *terms: [sum(k*w[i] for w, k in terms) for i in range(3)]
    zero = [ar.mpf(0)]*3
    if isinstance(s, Plane):
        jet = [comb((o, 1), (x, u), (y, v)), x, y, zero, zero, zero]
    elif isinstance(s, Cylinder):
        r = ar.mpf(s.radius)
        jet = [comb((o, 1), (e, r), (n, v)), comb((e1, r)), n, comb((e, -r)), zero, zero]
    elif isinstance(s, Cone):
        a = ar.mpf(s.half_angle)
        ca, sa = ar.cos(a), ar.sin(a)
        rho = ar.mpf(s.radius)+sa*v
        jet = [comb((o, 1), (e, rho), (n, ca*v)), comb((e1, rho)), comb((e, sa), (n, ca)), comb((e, -rho)),
               comb((e1, sa)), zero]
    elif isinstance(s, Sphere):
        r = ar.mpf(s.radius)
        cv, sv = ar.cos(v), ar.sin(v)
        jet = [comb((o, 1), (e, r*cv), (n, r*sv)), comb((e1, r*cv)), comb((e, -r*sv), (n, r*cv)),
               comb((e, -r*cv)), comb((e1, -r*sv)), comb((e, -r*cv), (n, -r*sv))]
    else:
        big, small = ar.mpf(s.major), ar.mpf(s.minor)
        cv, sv = ar.cos(v), ar.sin(v)
        rho = big+small*cv
        jet = [comb((o, 1), (e, rho), (n, small*sv)), comb((e1, rho)), comb((e, -small*sv), (n, small*cv)),
               comb((e, -rho)), comb((e1, -small*sv)), comb((e, -small*cv), (n, -small*sv))]
    return jet if interval else jet[:3]


def _upper(v):
    """An upper bound of the norm of a vector of intervals."""
    return mp.sqrt(sum(max(abs(x.a), abs(x.b))**2 for x in v))


def taylor_deviation(curve, s, p, forward, tol):
    """(low, high) of a use by second-order Taylor bounds on bisected pieces
    of the use's fraction, with interval jets: |D| <= |D(m)| + |D'(m)| h/2 +
    sup|D''| h^2/8. None on a periodic spline surface, as in the kernel."""
    if isinstance(s, BSplineSurface) and (s.u.periodic or s.v.periodic):
        return None
    cuts = {F(0), F(1)}
    if isinstance(p, BSpline2):
        cuts |= set(spline.knot_fractions(p))
    if isinstance(curve, BSpline3):
        cuts |= {k if forward else 1-k for k in spline.knot_fractions(curve)}
    cuts = sorted(cuts)
    sign = 1 if forward else -1
    edge_at = lambda t: t if forward else 1-t
    low, high = mp.mpf(0), mp.mpf(0)
    stack = [(t0, t1, 0) for t0, t1 in zip(cuts, cuts[1:])]
    while stack:
        t0, t1, depth = stack.pop()
        m, h = (t0+t1)/2, t1-t0
        cm, c1 = _point_jet(curve, mp_of(edge_at(m)))
        pm, p1 = _point_jet(p, mp_of(m))
        sm = _surface_jet(s, pm[0], pm[1], False)
        d0 = sub(cm, sm[0])
        d1 = [sign*a-b*p1[0]-c*p1[1] for a, b, c in zip(c1, sm[1], sm[2])]
        lo_edge, hi_edge = sorted((edge_at(t0), edge_at(t1)))
        box_t = mp.iv.mpf([mp_of(lo_edge), mp_of(hi_edge)])
        _, _, c2 = _interval_jet(curve, box_t)
        pv, pd1, pd2 = _interval_jet(p, mp.iv.mpf([mp_of(t0), mp_of(t1)]))
        jet = _surface_jet(s, pv[0], pv[1], True)
        if jet is None:
            return None
        _, su, sv, suu, suv, svv = jet
        d2 = [c-(a*pd1[0]**2+2*b*pd1[0]*pd1[1]+cc*pd1[1]**2+d*pd2[0]+e*pd2[1])
              for c, a, b, cc, d, e in zip(c2, suu, suv, svv, su, sv)]
        point = norm(d0)
        low = max(low, point)
        bound = point+norm(d1)*mp_of(h)/2+_upper(d2)*mp_of(h)**2/8
        if bound > mp.mpf(tol)/4 and depth < 14 and point <= mp.mpf(tol):
            m = (t0+t1)/2
            stack += [(t0, m, depth+1), (m, t1, depth+1)]
            continue
        high = max(high, bound)
    return low, high


def spline_face_flux(c, f):
    """The flux of a face on a spline surface by nested quadrature: G(u, v)
    = integral of S.(S_u x S_v) in v from the domain's start, then -loop
    integral of G du. None on a periodic surface, as in the kernel."""
    s = f.surface
    if s.u.periodic or s.v.periodic:
        return None
    _, _, (va, vb) = s.v.flat()
    vknots = [va]+[F(k) for k in s.v.knots if va < F(k) < vb]

    def flux_density(u, v):
        pnt, su, sv = spline.surface_jet_mp(s, u, v)
        return sum(a*b for a, b in zip(pnt, (su[1]*sv[2]-su[2]*sv[1], su[2]*sv[0]-su[0]*sv[2],
                                              su[0]*sv[1]-su[1]*sv[0])))

    def G(u, v):
        breaks = [mp_of(k) for k in vknots if mp_of(k) < v]+[v]
        return sum(gauss(lambda w: flux_density(u, w), a, b) for a, b in zip(breaks, breaks[1:]))
    total = mp.mpf(0)
    for lid in f.loops:
        loop = c.loops[lid]
        if loop.vertex is not None:
            continue
        pieces, chords = loop_points(c, loop)
        for k in loop.fins:
            p = c.fins[k].pcurve
            cuts = spline_cuts(None, s, p, True)
            for t0, t1 in zip(cuts, cuts[1:]):

                def integrand(t, p=p):
                    if isinstance(p, BSpline2):
                        q, d = spline.curve_jet_mp(p, t)
                        du = d[0]
                    else:
                        q = pcurve_point(p, t)
                        du = mp.mpf(p.end[0])-mp.mpf(p.start[0])
                    return -G(q[0], q[1])*du
                total += gauss(integrand, mp_of(t0), mp_of(t1))
        for a0, b0 in chords:
            if b0[0] != a0[0]:
                total += gauss(lambda t: -G(a0[0]+(b0[0]-a0[0])*t, a0[1]+(b0[1]-a0[1])*t)*(b0[0]-a0[0]), 0, 1)
    return total


def face_mass_terms(c, f, ref):
    """The fourteen mass integrals of a face relative to `ref`, in the
    kernel's order (volume, first moments, second moments, mixed moments
    xy, yz, zx, |N| and the face centre's three p|N|), by Green's theorem
    in UV with nested Gauss-Legendre quadrature of the exact surface jets:
    G(u, v) = integral of the integrands in v from the domain's start (0 on
    an analytic surface, broken at knots on a spline one), then -loop
    integral of G du. None where the kernel's route has no counterpart
    here: a periodic spline surface, a loop winding in v, a face whose
    windings do not cancel, an arc pcurve on a curved surface."""
    s = f.surface
    if isinstance(s, BSplineSurface):
        if s.u.periodic or s.v.periodic:
            return None
        _, _, (va, vb) = s.v.flat()
        breaks_v = [va]+[F(k) for k in s.v.knots if va < F(k) < vb]
    elif isinstance(s, (Plane, Cylinder, Cone)):
        breaks_v = [F(0)]
    else:
        return None
    loops = [c.loops[lid] for lid in f.loops if c.loops[lid].vertex is None]
    if any(lp.winding_v for lp in loops) or sum(lp.winding for lp in loops):
        return None
    ref = [mp.mpf(x) for x in ref]

    def density(u, v):
        pnt, su, sv = _surface_jet(s, u, v, False)
        p = [a-b for a, b in zip(pnt, ref)]
        n = [su[1]*sv[2]-su[2]*sv[1], su[2]*sv[0]-su[0]*sv[2], su[0]*sv[1]-su[1]*sv[0]]
        size = mp.sqrt(sum(x*x for x in n))
        sq = [x*x for x in p]
        return ([sum(a*b for a, b in zip(p, n))/3]
                + [sq[i]*n[i]/2 for i in range(3)]
                + [sq[i]*p[i]*n[i]/3 for i in range(3)]
                + [sq[i]*p[j]*n[i]/2 for i, j in ((0, 1), (1, 2), (2, 0))]
                + [size]+[x*size for x in p])

    def G(u, v):
        cuts = [mp_of(k) for k in breaks_v]
        # Below the first break the integral runs backwards.
        points = [x for x in cuts if x < v]+[v] if v >= cuts[0] else [cuts[0], v]
        total = [mp.mpf(0)]*14
        for a, b in zip(points, points[1:]):
            part = gauss_terms(lambda w: density(u, w), a, b)
            total = [x+y for x, y in zip(total, part)]
        return total

    total = [mp.mpf(0)]*14
    for loop in loops:
        pieces, chords = loop_points(c, loop)
        for k in loop.fins:
            p = c.fins[k].pcurve
            if isinstance(p, Arc2) and not isinstance(s, Plane):
                return None
            cuts = spline_cuts(None, s, p, True)
            for t0, t1 in zip(cuts, cuts[1:]):

                def integrand(t, p=p):
                    if isinstance(p, BSpline2):
                        q, d = spline.curve_jet_mp(p, t)
                        du = d[0]
                    elif isinstance(p, Line2):
                        q = pcurve_point(p, t)
                        du = mp.mpf(p.end[0])-mp.mpf(p.start[0])
                    else:
                        q = pcurve_point(p, t)
                        a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
                        du = -p.radius*mp.sin(a)*p.sweep
                    return [-x*du for x in G(q[0], q[1])]
                part = gauss_terms(integrand, mp_of(t0), mp_of(t1))
                total = [x+y for x, y in zip(total, part)]
        for a0, b0 in chords:
            if b0[0] != a0[0]:
                part = gauss_terms(lambda t: [-x*(b0[0]-a0[0]) for x in
                                              G(a0[0]+(b0[0]-a0[0])*t, a0[1]+(b0[1]-a0[1])*t)], 0, 1)
                total = [x+y for x, y in zip(total, part)]
    return total


def mass_properties(c):
    """Volume, surface area, centroid and the inertia matrix about it
    (products of inertia negated, as OCCT's MatrixOfInertia) of the solid
    regions, from `face_mass_terms`; None when a face is not integrated."""
    ref = c.vertices[0] if c.vertices else (0.0, 0.0, 0.0)
    region = [mp.mpf(0)]*10
    area = mp.mpf(0)
    counted = set()
    for r in c.regions:
        if r.kind != 'solid':
            continue
        for si in r.shells:
            for fi, side in c.shells[si].sides:
                terms = face_mass_terms(c, c.faces[fi], ref)
                if terms is None:
                    return None
                sign = 1 if side == 'F' else -1
                region = [x+sign*y for x, y in zip(region, terms[:10])]
                if fi not in counted:
                    counted.add(fi)
                    area += abs(terms[10])
    volume, mx, my, mz, sxx, syy, szz, sxy, syz, sxz = region
    centre = [mx/volume, my/volume, mz/volume]
    j = [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
    trace = sxx+syy+szz
    c2 = sum(x*x for x in centre)
    inertia = [[(trace if a == b else 0)-j[a][b]-volume*((c2 if a == b else 0)-centre[a]*centre[b])
                for b in range(3)] for a in range(3)]
    return {'volume': volume, 'area': area, 'centroid': [x+mp.mpf(r) for x, r in zip(centre, ref)],
            'inertia': inertia}


def gauss_terms(f, a, b):
    """`gauss` of a vector-valued integrand."""
    nodes = _gauss_nodes(mp.mp.prec)
    half, mid = (mp.mpf(b)-mp.mpf(a))/2, (mp.mpf(b)+mp.mpf(a))/2
    values = [f(mid+half*x) for x, _ in nodes]
    return [half*mp.fsum(w*v[k] for (_, w), v in zip(nodes, values)) for k in range(len(values[0]))]


def gauss(f, a, b):
    """Gauss-Legendre quadrature on [a, b] with mpmath's 24 fixed nodes:
    exact for the polynomial integrands of nonrational patches up to degree
    47."""
    nodes = _gauss_nodes(mp.mp.prec)
    half, mid = (mp.mpf(b)-mp.mpf(a))/2, (mp.mpf(b)+mp.mpf(a))/2
    return half*mp.fsum(w*f(mid+half*x) for x, w in nodes)


@functools.lru_cache(maxsize=None)
def _gauss_nodes(prec):
    from mpmath.calculus.quadrature import GaussLegendre
    return GaussLegendre(mp.mp).calc_nodes(4, prec)


def surface_distance(s, p):
    o, x, y, n = axes(s.frame)
    rel = sub(p, o)
    if isinstance(s, Plane):
        return abs(dot(rel, n))
    if isinstance(s, Cone):
        # Distance to the generatrix lines of both nappes in the meridian
        # half-plane: r c - R c - z s and r c + R c + z s.
        a = mp.mpf(s.half_angle)
        z = dot(rel, n)
        r = norm(sub(rel, mul(n, z)))
        R = mp.mpf(s.radius)
        return min(abs(r*mp.cos(a)-R*mp.cos(a)-z*mp.sin(a)), abs(r*mp.cos(a)+R*mp.cos(a)+z*mp.sin(a)))
    if isinstance(s, Sphere):
        return abs(norm(rel)-s.radius)
    if isinstance(s, Torus):
        z = dot(rel, n)
        rho = norm(sub(rel, mul(n, z)))
        return abs(mp.sqrt((rho-s.major)**2+z**2)-s.minor)
    radial = sub(rel, mul(n, dot(rel, n)))
    return abs(norm(radial)-s.radius)


def shell_point(c, si, vertex_ok):
    """A point of the shell: its first face's first fin's start."""
    f = c.faces[c.shells[si].sides[0][0]]
    if not f.loops:
        # A face without loops (a whole sphere): its point at (0, 0).
        return surface_point(f.surface, (mp.mpf(0), mp.mpf(0)))
    loop = c.loops[f.loops[0]]
    if loop.vertex is not None:
        return vec(c.vertices[loop.vertex])
    k = loop.fins[0]
    v = fin_vertices(c, k)[0]
    if v is not None:
        return vec(c.vertices[v]) if vertex_ok[v] else None
    return surface_point(f.surface, pcurve_point(c.fins[k].pcurve, 0))


def loop_points(c, loop):
    """Pcurves of an edge loop in order, and closure chords between them."""
    pieces = [c.fins[k].pcurve for k in loop.fins]
    chords = []
    for j, k in enumerate(loop.fins):
        a = pcurve_point(pieces[j], 1)
        b = pcurve_point(pieces[(j+1) % len(pieces)], 0)
        if j == len(pieces)-1:
            b = [b[0]+TAU*loop.winding, b[1]+TAU*loop.winding_v]
        chords.append((a, b))
    return pieces, chords


def loop_area(c, loop):
    """Signed area 1/2 * integral(u dv - v du), with closure chords."""
    from brep_reference import loop_area as seamed_area, Use
    uses = [Use(c.fins[k].edge, c.fins[k].forward, c.fins[k].pcurve) for k in loop.fins]
    return seamed_area(uses)


def periodic_area(c, loop):
    """-integral v du on the universal cover, closed by chords (the period
    shift included); seam segments would contribute nothing."""
    pieces, chords = loop_points(c, loop)
    total = mp.mpf(0)
    for p in pieces:
        def f(t, p=p):
            q = pcurve_point(p, t)
            if isinstance(p, BSpline2):
                du = spline.curve_jet_mp(p, t)[1][0]
            elif isinstance(p, Line2):
                du = mp.mpf(p.end[0])-mp.mpf(p.start[0])
            else:
                a = mp.mpf(p.start)+mp.mpf(p.sweep)*t
                du = -p.radius*mp.sin(a)*p.sweep
            return -q[1]*du
        total += mp.quad(f, [0, mp.mpf(1)/2, 1])
    for a, b in chords:
        total += -(a[1]+b[1])/2*(b[0]-a[0])
    return total


def periodic_area_v(c, loop):
    """Integral of u dv on the universal cover, closed by chords: the
    periodic area of a loop winding in v. Lines only on a torus."""
    pieces, chords = loop_points(c, loop)
    total = mp.mpf(0)
    for p in pieces:
        if isinstance(p, BSpline2):
            total += mp.quad(lambda t, p=p: (lambda q, d: q[0]*d[1])(*spline.curve_jet_mp(p, t)), [0, 1])
            continue
        total += (mp.mpf(p.start[0])+mp.mpf(p.end[0]))/2*(mp.mpf(p.end[1])-mp.mpf(p.start[1]))
    for a, b in chords:
        total += (a[0]+b[0])/2*(b[1]-a[1])
    return total


def winding(c, loop, point):
    from brep_reference import winding as seamed_winding, Use
    return seamed_winding([Use(c.fins[k].edge, c.fins[k].forward, c.fins[k].pcurve) for k in loop.fins], point)


def face_flux(c, f):
    """Integral over the face of S.(S_u x S_v) du dv = -loop integral of v f(u) du,
    f independent of v; loops closed by chords (their orientation carries the
    face sense)."""
    s = f.surface
    if isinstance(s, BSplineSurface):
        return spline_face_flux(c, f)
    # On an analytic curved surface the kernel integrates a spline pcurve
    # only on a cylinder and only at constant u (it adds nothing there).
    if not isinstance(s, Plane):
        for lid in f.loops:
            for k in c.loops[lid].fins:
                p = c.fins[k].pcurve
                if isinstance(p, BSpline2) and (not isinstance(s, Cylinder) or len({q[0] for q in p.poles}) > 1):
                    return None
    o, x, y, n = axes(s.frame)
    if isinstance(s, Cone):
        # S.(S_u x S_v) = rho(v) h(u); its v-antiderivative from the apex is
        # rho(v)^2 / (2 sin a) h(u), zero at the pole.
        a = mp.mpf(s.half_angle)
        sa, ca, R = mp.sin(a), mp.cos(a), mp.mpf(s.radius)
        ox, oy, on = dot(o, x), dot(o, y), dot(o, n)
        hu = lambda u: ca*R+ca*(ox*mp.cos(u)+oy*mp.sin(u))-sa*on
        G = lambda u, v: (R+sa*v)**2/(2*sa)*hu(u)
    elif isinstance(s, Sphere):
        # S.(S_u x S_v) = R^2 cos v (O.q + R), q the unit radial direction;
        # its v-antiderivative from the south pole, in closed form.
        R = mp.mpf(s.radius)
        ox, oy, on = dot(o, x), dot(o, y), dot(o, n)
        G = lambda u, v: R**2*((ox*mp.cos(u)+oy*mp.sin(u))*((v+mp.pi/2)/2+mp.sin(2*v)/4)
                               + on*(mp.sin(v)**2-1)/2+R*(mp.sin(v)+1))
    elif isinstance(s, Torus):
        # S.(S_u x S_v) = r (R + r cos v)(A(u) cos v + (O.n) sin v + R cos v + r)
        # with A(u) = O.x cos u + O.y sin u; its antiderivative in v from 0
        # (G) and in u from 0 (H), in closed form.
        R, r = mp.mpf(s.major), mp.mpf(s.minor)
        ox, oy, on = dot(o, x), dot(o, y), dot(o, n)
        A = lambda u: ox*mp.cos(u)+oy*mp.sin(u)
        c2 = lambda v: v/2+mp.sin(2*v)/4
        G = lambda u, v: r*(R*A(u)*mp.sin(v)+R*on*(1-mp.cos(v))+R*R*mp.sin(v)+R*r*v+r*A(u)*c2(v)
                            + r*on*mp.sin(v)**2/2+r*R*c2(v)+r*r*mp.sin(v))
        H = lambda u, v: r*((R*mp.cos(v)+r*mp.cos(v)**2)*(ox*mp.sin(u)+oy*(1-mp.cos(u)))
                            + u*(R*on*mp.sin(v)+R*R*mp.cos(v)+R*r+r*on*mp.cos(v)*mp.sin(v)
                                 + r*R*mp.cos(v)**2+r*r*mp.cos(v)))
        edge_loops = [c.loops[l] for l in f.loops if c.loops[l].vertex is None]
        if not edge_loops:
            # The whole torus: A integrates to nothing over a turn.
            whole = 6*mp.pi**2*R*r*r
            return whole if f.forward else -whole
        if any(l.winding_v for l in edge_loops):
            total = mp.mpf(0)
            for loop in edge_loops:
                pieces, chords = loop_points(c, loop)
                segments = [(p.start, p.end) for p in pieces]+chords
                for a0, b0 in segments:
                    a0, b0 = [mp.mpf(t) for t in a0], [mp.mpf(t) for t in b0]
                    if b0[1] != a0[1]:
                        total += mp.quad(lambda t: H(a0[0]+(b0[0]-a0[0])*t, a0[1]+(b0[1]-a0[1])*t)*(b0[1]-a0[1]),
                                         [0, mp.mpf(1)/2, 1])
            return total
    elif isinstance(s, Plane):
        h = dot(o, cross(x, y))
        G = lambda u, v, h=h: v*h
    else:
        r = mp.mpf(s.radius)
        a, b = dot(o, cross(x, n)), dot(o, cross(y, n))
        det = dot(x, cross(y, n))
        G = lambda u, v, r=r, a=a, b=b, det=det: v*(r*(-mp.sin(u)*a+mp.cos(u)*b)+r*r*det)
    total = mp.mpf(0)
    if isinstance(s, Sphere):
        # The north pole's line bounds the cover when the face reaches it:
        # a whole sphere, or a band closed there. From the south pole, G
        # vanishes at the south pole.
        north = lambda: mp.quad(lambda u: G(u, mp.pi/2), [0, mp.pi, 2*mp.pi])
        if not any(c.loops[l].vertex is None for l in f.loops):
            total += north() if f.forward else -north()
        elif pole_loop(c, f) is not None and pole_north(c, f):
            total += face_turns(c, f)*north()
    for lid in f.loops:
        loop = c.loops[lid]
        if loop.vertex is not None:
            continue
        pieces, chords = loop_points(c, loop)
        for p in pieces:
            def integrand(t, p=p):
                q = pcurve_point(p, t)
                if isinstance(p, BSpline2):
                    du = spline.curve_jet_mp(p, t)[1][0]
                elif isinstance(p, Line2):
                    du = mp.mpf(p.end[0])-mp.mpf(p.start[0])
                else:
                    ang = mp.mpf(p.start)+mp.mpf(p.sweep)*t
                    du = -p.radius*mp.sin(ang)*p.sweep
                return -G(q[0], q[1])*du
            total += mp.quad(integrand, [0, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, 1])
        for a0, b0 in chords:
            if b0[0] != a0[0]:
                total += mp.quad(lambda t: -G(a0[0]+(b0[0]-a0[0])*t, a0[1]+(b0[1]-a0[1])*t)*(b0[0]-a0[0]), [0, 1])
    return total


def two_sided(c, si):
    """The faces the shell lists on both sides (a sheet's, S6)."""
    sides = c.shells[si].sides
    return {f for f, side in sides if (f, opposite(side)) in sides}


def shell_flux(c, si):
    """None when a face's flux is not decided (S4d). A face the shell lists
    on both sides adds nothing (S6)."""
    total = mp.mpf(0)
    both = two_sided(c, si)
    for f, side in c.shells[si].sides:
        if f in both:
            continue
        flux = face_flux(c, c.faces[f])
        if flux is None:
            return None
        total += flux if side == 'F' else -flux
    return total


def inside(c, si, point):
    """Ray parity against the shell's faces, in high precision. On a cylinder a
    hit is inside the face when a +v ray from it in the universal cover
    crosses the face's loops an odd number of times."""
    direction = [mp.mpf('0.5773502691896257'), mp.mpf('0.6123724356957945'), mp.mpf('0.5400617248673217')]
    hits = 0
    for fi, _ in c.shells[si].sides:
        face = c.faces[fi]
        s = face.surface
        if isinstance(s, (Cone, Sphere, Torus)):
            return None
        o, x, y, n = axes(s.frame)
        rel = sub(point, o)
        if isinstance(s, Plane):
            den = dot(direction, n)
            if den != 0:
                t = -dot(rel, n)/den
                if t > 0:
                    q = add(rel, mul(direction, t))
                    uv = [dot(q, x), dot(q, y)]
                    loops = [c.loops[l] for l in face.loops if c.loops[l].vertex is None]
                    if sum(winding(c, l, uv) for l in loops):
                        hits += 1
            continue
        pr = sub(rel, mul(n, dot(rel, n)))
        pd = sub(direction, mul(n, dot(direction, n)))
        a, b, cc = dot(pd, pd), 2*dot(pr, pd), dot(pr, pr)-mp.mpf(s.radius)**2
        disc = b*b-4*a*cc
        if a == 0 or disc <= 0:
            continue
        for t in ((-b-mp.sqrt(disc))/(2*a), (-b+mp.sqrt(disc))/(2*a)):
            if t <= 0:
                continue
            q = add(rel, mul(direction, t))
            uv = [mp.atan2(dot(q, y), dot(q, x)), dot(q, n)]
            if cover_crossings(c, face, uv) % 2 == 1:
                hits += 1
    return hits % 2 == 1


def cover_crossings(c, face, uv):
    """Crossings of the +v ray from uv with every pcurve and closure chord of
    the face, taken with all u aliases k*TAU."""
    count = 0
    for lid in face.loops:
        loop = c.loops[lid]
        if loop.vertex is not None:
            continue
        pieces, chords = loop_points(c, loop)
        segments = [(pcurve_point(p, 0), pcurve_point(p, 1)) for p in pieces if isinstance(p, Line2)]
        segments += chords
        for a, b in segments:
            lo, hi = min(a[0], b[0]), max(a[0], b[0])
            for k in range(int(mp.floor((uv[0]-hi)/TAU))-1, int(mp.ceil((uv[0]-lo)/TAU))+2):
                u = uv[0]-k*TAU
                # Half-open in u: a crossing when exactly one end lies right of u.
                if (a[0] > u) == (b[0] > u):
                    continue
                v = a[1]+(b[1]-a[1])*(u-a[0])/(b[0]-a[0])
                if v > uv[1]:
                    count += 1
    return count
