//! S9c.1's results from the kept pieces: pieces of one input face kept the
//! same way joined into maximal faces across the edges they share (a full
//! circle's halves among them), edges joined where they run straight on
//! along one curve between the same faces (a seam's vertices dropped, a
//! closed curve left without vertices), each edge's curve rounded once
//! (lines, circular arcs, ellipse arcs on principal axes), pcurves on the
//! input faces' stored surfaces (exact projections where the surface's
//! frame does not hold the curve's axes), solids and their cavities, and
//! each slot's provenance by S9b.1's rules.
use super::super::polyhedra::Component;
use super::super::SlotPlan;
use super::graph::*;
use super::meet::Pos;
use super::model::*;
use super::num::*;
use crate::identity::{EntityId, EntityKind, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::split::rational_f64;
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Projection, Region, RegionId, RegionKind, Shell, ShellId, Side, Slot, Surface,
    TopologyParts, Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Point3, Result, Vec3};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{FRAC_PI_2, TAU};

/// A face loop's half-edges as result edges' uses.
type FinsOf<'a> = dyn Fn(&[(usize, bool)]) -> HLoop + 'a;

fn find(p: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while p[r] != r {
        r = p[r];
    }
    let mut y = x;
    while p[y] != r {
        let n = p[y];
        p[y] = r;
        y = n;
    }
    r
}

fn union(p: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(p, a), find(p, b));
    if ra != rb {
        p[ra.max(rb)] = ra.min(rb);
    }
}

/// A result face: its pieces, its input face (operand, a model face of
/// it), whether it keeps the input face's orientation, and its loops of
/// half-edges run about the result's outward normal.
struct RFace {
    pieces: Vec<usize>,
    op: usize,
    face: usize,
    behind: bool,
    loops: Vec<Vec<(usize, bool)>>,
}

/// A result edge: arrangement edges joined along one curve, each with its
/// direction along the chain, and its end vertices (none for a ring).
pub(super) struct REdge {
    pub(super) parts: Vec<(usize, bool)>,
    ends: Option<[usize; 2]>,
}

/// The names of input faces by where a result lies on them: a given
/// result's model face (S9e.1) holding several of its result faces names
/// each piece by the result face it lies in (`given_parts`).
struct Names<'a> {
    arr: &'a Arr,
    parts: BTreeMap<usize, EntityId>,
}

impl Names<'_> {
    fn multiple(&self, o: usize, f: usize) -> bool {
        self.arr.models[o]
            .given
            .as_ref()
            .is_some_and(|g| g.ids[f].len() > 1)
    }

    /// A piece's input face.
    fn piece(&self, pi: usize) -> EntityId {
        let p = &self.arr.pieces[pi];
        match self.parts.get(&pi) {
            Some(id) => *id,
            None => self.arr.models[p.op].faces[p.face].id,
        }
    }

    /// Input face `f` of operand `o` where an arrangement edge lies on it.
    fn by_edge(&self, o: usize, f: usize, gid: usize) -> EntityId {
        if self.multiple(o, f) {
            for (pi, p) in self.arr.pieces.iter().enumerate() {
                if p.op == o && p.face == f && p.loops.iter().flatten().any(|h| h.0 == gid) {
                    return self.piece(pi);
                }
            }
        }
        self.arr.models[o].faces[f].id
    }

    /// Input face `f` of operand `o` where an arrangement vertex lies on it.
    fn by_vertex(&self, o: usize, f: usize, v: usize) -> EntityId {
        if self.multiple(o, f) {
            for (pi, p) in self.arr.pieces.iter().enumerate() {
                if p.op == o
                    && p.face == f
                    && p.loops
                        .iter()
                        .flatten()
                        .any(|h| self.arr.edges[h.0].ends.contains(&v))
                {
                    return self.piece(pi);
                }
            }
        }
        self.arr.models[o].faces[f].id
    }
}

/// The result face each piece of a given result's model face holding
/// several result faces lies in: its pieces joined across the edges they
/// share (never across the model's edges between two result faces, which
/// no two pieces of one model face share), each group named by a model
/// edge on its boundary and the result face on its side.
fn given_parts(arr: &Arr) -> Result<BTreeMap<usize, EntityId>> {
    let mut out = BTreeMap::new();
    for o in 0..2 {
        let Some(g) = &arr.models[o].given else {
            continue;
        };
        for f in 0..g.ids.len() {
            if g.ids[f].len() < 2 {
                continue;
            }
            let ps: Vec<usize> = (0..arr.pieces.len())
                .filter(|&i| arr.pieces[i].op == o && arr.pieces[i].face == f)
                .collect();
            let mut parent: Vec<usize> = (0..ps.len()).collect();
            let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (k, &pi) in ps.iter().enumerate() {
                for &(gid, _) in arr.pieces[pi].loops.iter().flatten() {
                    users.entry(gid).or_default().push(k);
                }
            }
            for us in users.values() {
                for w in us.windows(2) {
                    union(&mut parent, w[0], w[1]);
                }
            }
            let mut named: BTreeMap<usize, EntityId> = BTreeMap::new();
            for (k, &pi) in ps.iter().enumerate() {
                for &(gid, _) in arr.pieces[pi].loops.iter().flatten() {
                    if let CurveRef::Edge(eo, ei) = arr.edges[gid].curve {
                        if let (true, Some(id)) = (eo == o, g.sides.get(&(ei, f))) {
                            let r = find(&mut parent, k);
                            named.entry(r).or_insert(*id);
                        }
                    }
                }
            }
            for (k, &pi) in ps.iter().enumerate() {
                let r = find(&mut parent, k);
                let id = named.get(&r).ok_or(Error::ComputationLimit(
                    "a given result's face part without its edges",
                ))?;
                out.insert(pi, *id);
            }
        }
    }
    Ok(out)
}

/// A result solid's provenance in its arrangement (S9e.1: a result given
/// to another Boolean is its arrangement's kept pieces): each face slot's
/// pieces, each edge slot's arrangement edges with their direction along
/// it, each vertex slot's arrangement vertex (none for a pole).
#[derive(Debug, Clone)]
pub(super) struct Made {
    pub(super) faces: Vec<Vec<usize>>,
    pub(super) edges: Vec<Vec<(usize, bool)>>,
    pub(super) vertices: Vec<Option<usize>>,
}

pub(super) fn assemble(arr: &Arr, op: Op2) -> Result<Vec<Component>> {
    Ok(assemble_made(arr, op)?.into_iter().map(|x| x.0).collect())
}

/// The result's solids with their provenance.
pub(super) fn assemble_made(arr: &Arr, op: Op2) -> Result<Vec<(Component, Made)>> {
    let kept: Vec<usize> = (0..arr.pieces.len())
        .filter(|&i| arr.pieces[i].keep)
        .collect();
    if kept.is_empty() {
        return Ok(Vec::new());
    }
    let names = Names {
        arr,
        parts: given_parts(arr)?,
    };
    // Pieces of one input face kept the same way, sharing an edge, join.
    let key = |i: usize| {
        let p = &arr.pieces[i];
        (p.op, arr.models[p.op].faces[p.face].id, p.behind)
    };
    // Pieces of A and B on one surface facing one way join too: their
    // result normals (each input face's, reversed unless `behind`) agree.
    let same_surface = |i: usize, j: usize| {
        let (p, q) = (&arr.pieces[i], &arr.pieces[j]);
        if p.op == q.op {
            return false;
        }
        let (a, b) = if p.op == 0 { (p, q) } else { (q, p) };
        if !arr.coinc.contains(&(a.face, b.face)) {
            return false;
        }
        let x = &arr.edges[a.loops[0][0].0].mid;
        let na = arr.models[0].normal_at(a.face, x);
        let nb = arr.models[1].normal_at(b.face, x);
        let facing = qqdot(&na, &nb).sign() == Ordering::Greater;
        facing == (a.behind == b.behind)
    };
    let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &i in &kept {
        for lp in &arr.pieces[i].loops {
            for &(g, _) in lp {
                users.entry(g).or_default().push(i);
            }
        }
    }
    let mut parent: Vec<usize> = (0..arr.pieces.len()).collect();
    for us in users.values() {
        for a in 0..us.len() {
            for b in a + 1..us.len() {
                if us[a] != us[b] && (key(us[a]) == key(us[b]) || same_surface(us[a], us[b])) {
                    union(&mut parent, us[a], us[b]);
                }
            }
        }
    }
    let mut group_of: BTreeMap<usize, usize> = BTreeMap::new();
    let mut faces: Vec<RFace> = Vec::new();
    for &i in &kept {
        let r = find(&mut parent, i);
        let gi = *group_of.entry(r).or_insert_with(|| {
            let p = &arr.pieces[i];
            faces.push(RFace {
                pieces: Vec::new(),
                op: p.op,
                face: p.face,
                behind: p.behind,
                loops: Vec::new(),
            });
            faces.len() - 1
        });
        faces[gi].pieces.push(i);
    }
    // Each face's loops: its pieces' half-edges less those it holds both
    // ways, linked again at their vertices. A piece of the other input on
    // the face's surface turned the other way (its input face's normal
    // against this one's: a reversed piece joining a kept one, S9e.2's
    // stack floors) runs its loops the other way, about the face's own
    // normal, and every vertex's way on is chosen about that normal.
    for rf in &mut faces {
        let mut face_of: BTreeMap<(usize, bool), usize> = BTreeMap::new();
        for &i in &rf.pieces {
            let same = arr.pieces[i].behind == rf.behind;
            for lp in &arr.pieces[i].loops {
                for &(g, d) in lp {
                    face_of.insert((g, if same { d } else { !d }), rf.face);
                }
            }
        }
        let hs: Vec<(usize, bool)> = face_of
            .keys()
            .copied()
            .filter(|&(g, d)| !face_of.contains_key(&(g, !d)))
            .collect();
        let loops = arr.relink(rf.op, &hs, &face_of)?;
        rf.loops = if rf.behind {
            loops
        } else {
            loops
                .into_iter()
                .map(|l| l.into_iter().rev().map(|(g, d)| (g, !d)).collect())
                .collect()
        };
    }
    // Every edge used once each way.
    let mut uses: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
    for (fi, rf) in faces.iter().enumerate() {
        for lp in &rf.loops {
            for &(g, d) in lp {
                uses.entry(g).or_default().push((fi, d));
            }
        }
    }
    for us in uses.values() {
        let fwd = us.iter().filter(|u| u.1).count();
        match (us.len(), fwd) {
            (2, 1) => {}
            (4, 2) => return Err(Error::Degenerate("a result touching itself along an edge")),
            _ => {
                return Err(Error::InvalidTopology(
                    "an open Boolean of arcs in any position",
                ))
            }
        }
    }
    // Vertices where two edges of one curve meet between the same faces.
    let mut incident: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for &g in uses.keys() {
        for v in arr.edges[g].ends {
            incident.entry(v).or_default().insert(g);
        }
    }
    let faces_of = |g: usize| -> BTreeSet<usize> { uses[&g].iter().map(|u| u.0).collect() };
    let removable = |v: usize| -> bool {
        // A pole stays: its sections' pcurves turn there.
        if matches!(arr.vx[v].key, VKey::Pole(..)) {
            return false;
        }
        let gs: Vec<usize> = incident[&v].iter().copied().collect();
        match gs.len() {
            2 => {
                same_curve(&arr.edges[gs[0]].crv, &arr.edges[gs[1]].crv)
                    && faces_of(gs[0]) == faces_of(gs[1])
            }
            // A closed curve's piece from and back to this vertex alone.
            1 => arr.edges[gs[0]].ends[0] == arr.edges[gs[0]].ends[1],
            _ => false,
        }
    };
    let removed: BTreeSet<usize> = incident.keys().copied().filter(|&v| removable(v)).collect();
    // Chains through removed vertices.
    let mut redges: Vec<REdge> = Vec::new();
    let mut chain_of: BTreeMap<usize, (usize, bool)> = BTreeMap::new();
    let other_at =
        |v: usize, g: usize| -> usize { *incident[&v].iter().find(|&&x| x != g).unwrap_or(&g) };
    for &g0 in uses.keys() {
        if chain_of.contains_key(&g0) {
            continue;
        }
        // Back to the chain's start.
        let (mut g, mut d) = (g0, true);
        let mut steps = 0;
        loop {
            let v = if d {
                arr.edges[g].ends[0]
            } else {
                arr.edges[g].ends[1]
            };
            if !removed.contains(&v) || steps > uses.len() {
                break;
            }
            let h = other_at(v, g);
            let hd = arr.edges[h].ends[1] == v;
            if (h, hd) == (g0, true) {
                break;
            }
            g = h;
            d = hd;
            steps += 1;
        }
        let first = (g, d);
        let mut parts = vec![first];
        let ring;
        loop {
            let (g, d) = *parts.last().expect("a part");
            let v = if d {
                arr.edges[g].ends[1]
            } else {
                arr.edges[g].ends[0]
            };
            if !removed.contains(&v) {
                ring = false;
                break;
            }
            let h = other_at(v, g);
            let hd = arr.edges[h].ends[0] == v;
            if (h, hd) == first {
                ring = true;
                break;
            }
            parts.push((h, hd));
        }
        let ends = if ring {
            None
        } else {
            let s = if first.1 {
                arr.edges[first.0].ends[0]
            } else {
                arr.edges[first.0].ends[1]
            };
            let (lg, ld) = *parts.last().expect("a part");
            let e = if ld {
                arr.edges[lg].ends[1]
            } else {
                arr.edges[lg].ends[0]
            };
            Some([s, e])
        };
        let ri = redges.len();
        for &(g, d) in &parts {
            chain_of.insert(g, (ri, d));
        }
        redges.push(REdge { parts, ends });
    }
    // Faces' loops as result edges' uses.
    let fins_of = |lp: &[(usize, bool)]| -> Vec<(usize, bool)> {
        let mut out: Vec<(usize, bool)> = Vec::new();
        for &(g, d) in lp {
            let (ri, cd) = chain_of[&g];
            let use_ = (ri, cd == d);
            if out.last() != Some(&use_) {
                out.push(use_);
            }
        }
        while out.len() > 1 && out.first() == out.last() {
            out.pop();
        }
        out
    };
    // Shells: faces joined by edges.
    let mut sp: Vec<usize> = (0..faces.len()).collect();
    for us in uses.values() {
        union(&mut sp, us[0].0, us[1].0);
    }
    let mut shells: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for fi in 0..faces.len() {
        let r = find(&mut sp, fi);
        shells.entry(r).or_default().push(fi);
    }
    let shells: Vec<Vec<usize>> = shells.into_values().collect();
    // Shells meeting at a vertex touch.
    let mut shell_of_vertex: BTreeMap<usize, usize> = BTreeMap::new();
    for (si, sh) in shells.iter().enumerate() {
        for &fi in sh {
            for lp in &faces[fi].loops {
                for &(g, _) in lp {
                    for v in arr.edges[g].ends {
                        if removed.contains(&v) {
                            continue;
                        }
                        if let Some(&s) = shell_of_vertex.get(&v) {
                            if s != si {
                                return Err(Error::Degenerate("solids touching at a vertex"));
                            }
                        }
                        shell_of_vertex.insert(v, si);
                    }
                }
            }
        }
    }
    // Every result vertex apart from the others by the resolution.
    let tol = arr.models[0].tolerance.linear();
    let points: BTreeMap<usize, Point3> = shell_of_vertex
        .keys()
        .map(|&v| {
            let p = qv_f64(&arr.vx[v].p);
            (v, Point3::new(p[0], p[1], p[2]))
        })
        .collect();
    let pv: Vec<(&usize, &Point3)> = points.iter().collect();
    for i in 0..pv.len() {
        for j in i + 1..pv.len() {
            if (*pv[i].1 - *pv[j].1).length() <= tol {
                return Err(Error::Degenerate("a result thinner than the resolution"));
            }
        }
    }
    let info: BTreeMap<EntityId, (Operand, Role)> = arr
        .models
        .iter()
        .flat_map(|m| m.info.iter().map(|(k, v)| (*k, *v)))
        .collect();
    // A shell of the tool's faces alone in a cut is a cavity; so is any
    // other shell bounding a void (S9d.4c: a band's end disc and inner wall
    // under another input's face over its hole, in a fuse), which the
    // validator's certified flux finds turned inward when the shell is
    // built alone. A void bounded by one input's faces alone would be a
    // cavity of that input: only shells of both inputs' faces are tried,
    // and only beside another shell.
    let mut is_cavity: Vec<bool> = shells
        .iter()
        .map(|sh| op == Op2::Cut && sh.iter().all(|&fi| faces[fi].op == 1))
        .collect();
    if shells.len() > 1 {
        for si in 0..shells.len() {
            let both = shells[si].iter().any(|&fi| faces[fi].op == 0)
                && shells[si].iter().any(|&fi| faces[fi].op == 1);
            if is_cavity[si] || !both {
                continue;
            }
            let (alone, _) = build_component(
                arr,
                &names,
                op,
                &faces,
                &redges,
                &fins_of,
                &points,
                &shells[si],
                &BTreeSet::new(),
                &info,
            )?;
            if let Err(issues) = crate::topology::Topology::from_parts(
                alone.parts.with_measured_enclosures(),
                arr.models[0].tolerance,
            ) {
                if issues
                    .iter()
                    .any(|i| i.kind == crate::topology::IssueKind::ShellOrientation)
                {
                    is_cavity[si] = true;
                }
            }
        }
    }
    let outers: Vec<usize> = (0..shells.len()).filter(|&s| !is_cavity[s]).collect();
    let cavities: Vec<usize> = (0..shells.len()).filter(|&s| is_cavity[s]).collect();
    if !cavities.is_empty() && outers.len() != 1 {
        return Err(Error::OutOfDomain("a cavity among several solids (S9c)"));
    }
    let mut out = Vec::new();
    for &si in &outers {
        let mut all: Vec<usize> = shells[si].clone();
        let mut inner: BTreeSet<usize> = BTreeSet::new();
        for &c in &cavities {
            all.extend(&shells[c]);
            inner.extend(&shells[c]);
        }
        out.push(build_component(
            arr, &names, op, &faces, &redges, &fins_of, &points, &all, &inner, &info,
        )?);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn build_component(
    arr: &Arr,
    names: &Names<'_>,
    op: Op2,
    faces: &[RFace],
    redges: &[REdge],
    fins_of: &FinsOf<'_>,
    points: &BTreeMap<usize, Point3>,
    all: &[usize],
    inner: &BTreeSet<usize>,
    info: &BTreeMap<EntityId, (Operand, Role)>,
) -> Result<(Component, Made)> {
    let mut p = TopologyParts::default();
    // Edges and vertices used.
    let mut used_edges: BTreeSet<usize> = BTreeSet::new();
    let face_fins: Vec<Vec<Vec<(usize, bool)>>> = all
        .iter()
        .map(|&fi| faces[fi].loops.iter().map(|l| fins_of(l)).collect())
        .collect();
    for loops in &face_fins {
        for l in loops {
            used_edges.extend(l.iter().map(|u| u.0));
        }
    }
    let mut verts: BTreeSet<usize> = BTreeSet::new();
    for &ri in &used_edges {
        if let Some(e) = redges[ri].ends {
            verts.extend(e);
        }
    }
    let vertex_id: BTreeMap<usize, VertexId> = verts
        .iter()
        .enumerate()
        .map(|(i, v)| (*v, VertexId(i)))
        .collect();
    for v in &verts {
        p.vertices.push(Vertex {
            position: points[v],
            enclosure: None,
        });
    }
    let edge_id: BTreeMap<usize, EdgeId> = used_edges
        .iter()
        .enumerate()
        .map(|(i, e)| (*e, EdgeId(i)))
        .collect();
    for &ri in &used_edges {
        let e = &redges[ri];
        p.edges.push(Edge {
            start: e.ends.map(|x| vertex_id[&x[0]]),
            end: e.ends.map(|x| vertex_id[&x[1]]),
            curve: curve3(arr, e, points)?,
            fins: Vec::new(),
        });
    }
    let has_cavity = !inner.is_empty();
    // Poles added to sphere faces whose loops wind once: (vertex, operand,
    // model face).
    let mut poles: Vec<(VertexId, usize, usize)> = Vec::new();
    for (k, &fi) in all.iter().enumerate() {
        let rf = &faces[fi];
        let mface = &arr.models[rf.op].faces[rf.face];
        let surface = mface.stored.clone();
        let sense = if rf.behind {
            mface.sense
        } else {
            flip(mface.sense)
        };
        let mut loop_ids = Vec::new();
        for l in &face_fins[k] {
            let fins = loop_fins(&p, l, &edge_id, &surface)?;
            let mut fids = Vec::new();
            for fin in fins {
                let id = FinId(p.fins.len());
                p.edges[fin.edge.0].fins.push(id);
                p.fins.push(fin);
                fids.push(id);
            }
            let winding = match &surface {
                Surface::Cylinder { .. } | Surface::Sphere { .. } | Surface::Cone { .. } => {
                    [turns(&p, &fids, &surface), 0]
                }
                // A torus's loops wind in either angle (S9d.4a).
                Surface::Torus { .. } => [turns(&p, &fids, &surface), turns_v(&p, &fids)],
                _ => [0, 0],
            };
            // The loop's area in the surface's parameters (its pcurves).
            let mut pts: Vec<Point2> = Vec::new();
            for f in &fids {
                let pc = &p.fins[f.0].pcurve;
                for i in 0..16 {
                    pts.push(pc.point(i as f64 / 16.0));
                }
            }
            let area = (0..pts.len())
                .map(|i| {
                    let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                    a.x * b.y - a.y * b.x
                })
                .sum::<f64>()
                .abs();
            loop_ids.push((LoopId(p.loops.len()), winding[0] != 0, area));
            p.loops.push(Loop::Edges {
                fins: fids,
                winding,
            });
        }
        // The outer loop first: loops round the axis, then by area.
        loop_ids.sort_by(|a, b| b.1.cmp(&a.1).then(b.2.total_cmp(&a.2)));
        let loop_ids: Vec<LoopId> = loop_ids.into_iter().map(|x| x.0).collect();
        // On a cylinder, each other loop lifted by whole turns to lie with
        // the first (each pcurve's lift starts from its point's angle).
        if matches!(
            surface,
            Surface::Cylinder { .. }
                | Surface::Sphere { .. }
                | Surface::Cone { .. }
                | Surface::Torus { .. }
        ) && loop_ids.len() > 1
        {
            let mean_u = |p: &TopologyParts, l: LoopId| -> f64 {
                let Loop::Edges { fins, .. } = &p.loops[l.0] else {
                    return 0.0;
                };
                let us: Vec<f64> = fins
                    .iter()
                    .flat_map(|f| (0..8).map(move |i| (f, i as f64 / 8.0)))
                    .map(|(f, t)| p.fins[f.0].pcurve.point(t).x)
                    .collect();
                us.iter().sum::<f64>() / us.len() as f64
            };
            // A torus's loops (S9d.4a): a reference winding loop fixes the
            // sheet, every other loop placed on the side its material lies
            // (above a loop running +u in v, below one running -u; below a
            // loop running +v in u, above one running -v), within a turn.
            if matches!(surface, Surface::Torus { .. }) {
                torus_sheets(&mut p, &loop_ids, sense == Orientation::Reversed);
                torus_holes_in(&mut p, &loop_ids);
            }
            let target = mean_u(&p, loop_ids[0]);
            let torus = matches!(surface, Surface::Torus { .. });
            for &l in loop_ids[1..].iter().filter(|_| !torus) {
                let turns = ((target - mean_u(&p, l)) / TAU).round();
                if turns == 0.0 {
                    continue;
                }
                let Loop::Edges { fins, .. } = &p.loops[l.0] else {
                    continue;
                };
                for f in fins.clone() {
                    let k = TAU * turns;
                    match &mut p.fins[f.0].pcurve {
                        Curve2::Projection(pr) => {
                            for lift in &mut pr.lifts {
                                lift.x += k;
                            }
                        }
                        Curve2::LineSegment { start, end } => {
                            start.x += k;
                            end.x += k;
                        }
                        Curve2::Sinusoid { start, .. } => *start += k,
                        _ => {}
                    }
                }
            }
        }
        // A sphere's face whose loops wind once in all closes at a pole, a
        // vertex loop (S3's caps): a band winding `+u` on a forward face
        // closes at the north pole.
        let mut loop_ids = loop_ids;
        if let Surface::Sphere { frame, radius } = &surface {
            let total: i32 = loop_ids
                .iter()
                .map(|l| match &p.loops[l.0] {
                    Loop::Edges { winding, .. } => winding[0],
                    Loop::Vertex(_) => 0,
                })
                .sum();
            if total.abs() == 1 {
                let north = (total == 1) == (sense == Orientation::Forward);
                let at = frame.point(Point2::default(), if north { *radius } else { -*radius });
                let vid = VertexId(p.vertices.len());
                p.vertices.push(Vertex {
                    position: at,
                    enclosure: None,
                });
                p.loops.push(Loop::Vertex(vid));
                loop_ids.push(LoopId(p.loops.len() - 1));
                poles.push((vid, rf.op, rf.face));
            }
        }
        // A cone's wall whose loops wind round its axis closes at its apex,
        // a vertex loop (S3's cones; S9d.3a).
        if let Surface::Cone { .. } = &surface {
            let total: i32 = loop_ids
                .iter()
                .map(|l| match &p.loops[l.0] {
                    Loop::Edges { winding, .. } => winding[0],
                    Loop::Vertex(_) => 0,
                })
                .sum();
            if total != 0 {
                let (m, _) = arr.models[rf.op].view(rf.face);
                let apex = m
                    .funnel
                    .as_ref()
                    .and_then(|f| f.apex)
                    .filter(|_| total.abs() == 1)
                    .ok_or(Error::InvalidTopology(
                        "a cone's wall winding without an apex",
                    ))?;
                let at = qv_f64(&m.verts[apex].p);
                let vid = VertexId(p.vertices.len());
                p.vertices.push(Vertex {
                    position: Point3::new(at[0], at[1], at[2]),
                    enclosure: None,
                });
                p.loops.push(Loop::Vertex(vid));
                loop_ids.push(LoopId(p.loops.len() - 1));
                poles.push((vid, rf.op, rf.face));
            }
        }
        let is_inner = inner.contains(&fi);
        p.faces.push(Face {
            surface,
            sense,
            loops: loop_ids,
            front: if is_inner { ShellId(2) } else { ShellId(0) },
            back: if is_inner { ShellId(3) } else { ShellId(1) },
            enclosure: None,
        });
    }
    let sides = |want_inner: bool, side: Side| -> Vec<(FaceId, Side)> {
        all.iter()
            .enumerate()
            .filter(|(_, fi)| inner.contains(fi) == want_inner)
            .map(|(k, _)| (FaceId(k), side))
            .collect()
    };
    let shell = |region: usize, sides: Vec<(FaceId, Side)>| Shell {
        region: RegionId(region),
        sides,
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    };
    p.shells = vec![
        shell(1, sides(false, Side::Front)),
        shell(0, sides(false, Side::Back)),
    ];
    p.regions = vec![
        Region {
            kind: RegionKind::Void,
            shells: vec![ShellId(1)],
        },
        Region {
            kind: RegionKind::Solid,
            shells: vec![ShellId(0)],
        },
    ];
    if has_cavity {
        p.shells.push(shell(1, sides(true, Side::Front)));
        p.shells.push(shell(2, sides(true, Side::Back)));
        p.regions[1].shells.push(ShellId(2));
        p.regions.push(Region {
            kind: RegionKind::Void,
            shells: vec![ShellId(3)],
        });
    }
    // Plans.
    let tool = |o: usize| op == Op2::Cut && o == 1;
    let tidy = |mut c: Vec<EntityId>, mut t: Vec<EntityId>| {
        c.sort();
        c.dedup();
        t.sort();
        t.dedup();
        t.retain(|x| !c.contains(x));
        (c, t)
    };
    let operand = |id: &EntityId| info.get(id).map(|i| i.0);
    let role_of =
        |ids: &[EntityId], new: Role| ids.first().and_then(|id| info.get(id)).map_or(new, |i| i.1);
    let mut plans: Vec<SlotPlan> = Vec::new();
    let mut members: BTreeSet<Operand> = BTreeSet::new();
    for (k, &fi) in all.iter().enumerate() {
        let rf = &faces[fi];
        let (mut c, mut t) = (Vec::new(), Vec::new());
        for &pi in &rf.pieces {
            let piece = &arr.pieces[pi];
            let id = names.piece(pi);
            if piece.behind && !tool(piece.op) {
                c.push(id);
            } else {
                t.push(id);
            }
        }
        members.extend(c.iter().filter_map(operand));
        let (c, t) = tidy(c, t);
        let role = role_of(&c, Role::CutFace);
        plans.push((Slot::Face(FaceId(k)), c, t, EntityKind::Face, role));
    }
    for (&ri, &eid) in &edge_id {
        let (mut c, mut t) = (Vec::new(), Vec::new());
        for &(g, _) in &redges[ri].parts {
            match arr.edges[g].curve {
                CurveRef::Edge(o, ei) => {
                    let me = &arr.models[o].edges[ei];
                    match me.id {
                        Some(id) if tool(o) => t.push(id),
                        Some(id) => c.push(id),
                        None => t.extend(me.faces.iter().map(|&f| names.by_edge(o, f, g))),
                    }
                }
                CurveRef::Section(si, _) => {
                    let s = &arr.secs[si];
                    t.push(names.by_edge(0, s.fa, g));
                    t.push(names.by_edge(1, s.fb, g));
                }
            }
        }
        let (c, t) = tidy(c, t);
        let role = role_of(&c, Role::CutEdge);
        plans.push((Slot::Edge(eid), c, t, EntityKind::Edge, role));
    }
    // A pole continues the input's pole vertex where it had one, else it is
    // generated from the sphere's face.
    for &(vid, o, f) in &poles {
        let m = &arr.models[o];
        let at = p.vertices[vid.0].position.to_array();
        let tol = m.tolerance.linear();
        let input = m.verts.iter().find_map(|v| {
            let q = qv_f64(&v.p);
            ((0..3).all(|k| (q[k] - at[k]).abs() <= tol))
                .then_some(v.id)
                .flatten()
        });
        let (c, t) = match input {
            Some(id) if tool(o) => (Vec::new(), vec![id]),
            Some(id) => (vec![id], Vec::new()),
            None => (Vec::new(), vec![m.faces[f].id]),
        };
        let (c, t) = tidy(c, t);
        let role = role_of(&c, Role::CutVertex);
        plans.push((Slot::Vertex(vid), c, t, EntityKind::Vertex, role));
    }
    for (&v, &vid) in &vertex_id {
        let (mut c, mut t) = (Vec::new(), Vec::new());
        match &arr.vx[v].key {
            VKey::Input(o, i) => match arr.models[*o].verts[*i].id {
                Some(id) if tool(*o) => t.push(id),
                Some(id) => c.push(id),
                None => {
                    // A given result's vertex inside one of its edges (a
                    // ring's seam, S9e.1) lies on that edge.
                    let m = &arr.models[*o];
                    let on: Vec<EntityId> = match &m.given {
                        Some(_) => m
                            .edges
                            .iter()
                            .filter(|e| e.start == *i || e.end == *i)
                            .filter_map(|e| e.id)
                            .collect(),
                        None => Vec::new(),
                    };
                    if on.is_empty() {
                        for &(o2, f) in &arr.vx[v].faces {
                            t.push(names.by_vertex(o2, f, v));
                        }
                    } else {
                        t.extend(on);
                    }
                }
            },
            VKey::Pierce(o, ei, g, _) => {
                let me = &arr.models[*o].edges[*ei];
                match me.id {
                    Some(id) => t.push(id),
                    None => t.extend(me.faces.iter().map(|&f| names.by_vertex(*o, f, v))),
                }
                t.push(names.by_vertex(1 - o, *g, v));
            }
            VKey::Cross(fa, fb, _) => {
                t.push(names.by_vertex(0, *fa, v));
                t.push(names.by_vertex(1, *fb, v));
            }
            VKey::Ring(si, _) | VKey::Pole(si, _) => {
                let s = &arr.secs[*si];
                t.push(names.by_vertex(0, s.fa, v));
                t.push(names.by_vertex(1, s.fb, v));
            }
        }
        let (c, t) = tidy(c, t);
        let role = role_of(&c, Role::CutVertex);
        plans.push((Slot::Vertex(vid), c, t, EntityKind::Vertex, role));
    }
    if op == Op2::Cut {
        members = BTreeSet::from([Operand::A]);
    }
    plans.push((
        Slot::Region(RegionId(1)),
        members
            .iter()
            .map(|o| arr.models[usize::from(*o == Operand::B)].region)
            .collect(),
        Vec::new(),
        EntityKind::Region,
        Role::Region,
    ));
    if has_cavity {
        plans.push((
            Slot::Region(RegionId(2)),
            Vec::new(),
            vec![arr.models[1].region],
            EntityKind::Region,
            Role::Region,
        ));
    }
    let mut vertices: Vec<Option<usize>> = vec![None; p.vertices.len()];
    for (&v, &vid) in &vertex_id {
        vertices[vid.0] = Some(v);
    }
    let made = Made {
        faces: all.iter().map(|&fi| faces[fi].pieces.clone()).collect(),
        edges: used_edges
            .iter()
            .map(|&ri| redges[ri].parts.clone())
            .collect(),
        vertices,
    };
    Ok((Component { parts: p, plans }, made))
}

/// Whether two curves are one (exactly): a full circle's halves, or a
/// plane's sections of a cylinder's halves, run on as one edge.
fn same_curve(a: &Crv, b: &Crv) -> bool {
    match (a, b) {
        (
            Crv::Conic { c, a: x, b: y },
            Crv::Conic {
                c: c2,
                a: x2,
                b: y2,
            },
        ) => c == c2 && x == x2 && y == y2,
        (Crv::Line { p, d }, Crv::Line { p: p2, d: d2 }) => {
            is_zero(&cross(d, d2)) && super::graph::on_line(p, d, p2)
        }
        (Crv::Meet(x), Crv::Meet(y)) => x == y,
        (Crv::Circle(x), Crv::Circle(y)) => x == y,
        (Crv::Rise(x), Crv::Rise(y)) => x == y,
        (Crv::Cone(x), Crv::Cone(y)) => x == y,
        (Crv::Torus(x), Crv::Torus(y)) => x == y,
        (Crv::Toric(x), Crv::Toric(y)) => x == y,
        _ => false,
    }
}

fn flip(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
    }
}

/// The binary64 angle of a place on a conic.
fn angle_of(pos: &Pos) -> f64 {
    let Pos::Ang(cs) = pos else {
        unreachable!("a conic's place")
    };
    cs[1].to_f64().atan2(cs[0].to_f64())
}

/// A result edge's curve, rounded once.
fn curve3(arr: &Arr, e: &REdge, points: &BTreeMap<usize, Point3>) -> Result<Curve3> {
    let (g0, d0) = e.parts[0];
    let (gl, dl) = *e.parts.last().expect("a part");
    let first = &arr.edges[g0];
    let last = &arr.edges[gl];
    // A torus segment's or wedge's rim (S9d.4b.1): on the input's stored
    // circle, running with its angle (the model's `u` or `v`).
    if let CurveRef::Edge(o, ei) = first.curve {
        let m = &arr.models[o];
        if let (EdgeKind::Rim(high, _), Some(ring)) = (m.edges[ei].kind, &m.ring) {
            let Some(Curve3::Circle { frame, radius }) = &ring.rims[usize::from(high)] else {
                return Err(Error::InvalidTopology("a torus rim off a circle"));
            };
            let Some([s, t]) = e.ends else {
                return Ok(Curve3::Circle {
                    frame: *frame,
                    radius: *radius,
                });
            };
            let angle = |p: Point3| {
                let [x, y, _] = frame.coordinates(p);
                y.atan2(x)
            };
            let (t0, t1) = (angle(points[&s]), angle(points[&t]));
            let with = first.with == d0;
            let mut sweep = if with { t1 - t0 } else { t0 - t1 }.rem_euclid(TAU);
            if sweep == 0.0 {
                sweep = TAU;
            }
            return Ok(Curve3::CircularArc {
                frame: *frame,
                radius: *radius,
                start_angle: t0,
                sweep_angle: if with { sweep } else { -sweep },
            });
        }
    }
    // S9e.3a: a result edge over the whole of a given edge of a procedural
    // curve keeps that edge's stored curve; a piece of a meeting of two
    // curved faces is S9e.3b's (its meetings are refused before).
    if let CurveRef::Edge(o, ei) = first.curve {
        if let (Some(g), true) = (&arr.models[o].given, super::chain::procedural(&first.crv)) {
            let whole = super::chain::whole(arr, e, o, ei);
            let ends = e.ends.map(|[s, t]| (points[&s], points[&t]));
            if let Some(c) = g.curves[ei]
                .as_ref()
                .filter(|_| whole)
                .and_then(|c| super::chain::stored(c, ends))
            {
                return Ok(c);
            }
            if !matches!(first.crv, Crv::Cone(_) | Crv::Torus(_)) {
                return Err(Error::OutOfDomain(super::chain::S9E3B));
            }
        }
    }
    match &first.crv {
        Crv::Rise(m) => {
            // Over the carrier's stored cylinder, heights from its origin.
            let CurveRef::Section(si, _) = first.curve else {
                unreachable!("a rise is a section")
            };
            let s = &arr.secs[si];
            let stored = [
                &arr.models[0].faces[s.fa].stored,
                &arr.models[1].faces[s.fb].stored,
            ];
            // Over a cone's stored surface too (S9d.3b.2).
            let (frame, radius, half_angle) = match stored[m.carrier] {
                Surface::Cylinder { frame, radius } => (frame, radius, 0.0),
                Surface::Cone {
                    frame,
                    radius,
                    half_angle,
                } => (frame, radius, *half_angle),
                _ => return Err(Error::InvalidTopology("a rise off a ruled surface")),
            };
            let (Pos::T(w0), Pos::T(w1)) = (
                if d0 { &first.pos[0] } else { &first.pos[1] },
                if dl { &last.pos[1] } else { &last.pos[0] },
            ) else {
                unreachable!("a rise's places")
            };
            let fl = |x: &V| x.clone().map(|y| rational_f64(&y));
            // The carrier's primitive model (a given face's view, S9e.3a).
            let n = super::chain::curve_model(arr, first.curve, m.carrier, |_| true)
                .f
                .n
                .clone();
            // The stored frame's origin's height along the model's axis,
            // from the model's circle centre at height 0.
            let origin = frame.origin().to_array();
            let nf = fl(&n);
            let cf = fl(&m.o_model());
            let offset: f64 = (0..3).map(|j| (origin[j] - cf[j]) * nf[j]).sum();
            let centre = fl(&m.c);
            Ok(Curve3::Rise(Box::new(crate::topology::Rise {
                frame: *frame,
                radius: *radius,
                half_angle,
                centre: Point3::new(centre[0], centre[1], centre[2]),
                sphere_radius: rational_f64(&m.rr),
                sign: if m.plus { 1.0 } else { -1.0 },
                start: w0.to_f64() - offset,
                sweep: w1.to_f64() - w0.to_f64(),
            })))
        }
        Crv::Cone(c) => cone_curve3(
            arr,
            super::chain::curve_model(arr, first.curve, c.carrier, |v| v.funnel.is_some()),
            e,
            c,
            points,
        ),
        Crv::Toric(c) => {
            // S9d.4b.2: on the torus's and the quadric's (or the other
            // torus's, S9d.4b.2b) stored surfaces.
            let CurveRef::Section(si, _) = first.curve else {
                unreachable!("a torus meeting is a section")
            };
            let s = &arr.secs[si];
            let stored = [
                &arr.models[0].faces[s.fa].stored,
                &arr.models[1].faces[s.fb].stored,
            ];
            let Surface::Torus {
                frame,
                major,
                minor,
            } = stored[c.carrier]
            else {
                return Err(Error::InvalidTopology("a torus meeting off a torus"));
            };
            let (other, other_radius, other_sphere, other_half_angle, other_minor) =
                match stored[1 - c.carrier] {
                    Surface::Cylinder { frame, radius } => (*frame, *radius, false, 0.0, 0.0),
                    Surface::Sphere { frame, radius } => (*frame, *radius, true, 0.0, 0.0),
                    Surface::Cone {
                        frame,
                        radius,
                        half_angle,
                    } => (*frame, *radius, false, *half_angle, 0.0),
                    Surface::Torus {
                        frame,
                        major,
                        minor,
                    } => (*frame, *major, false, 0.0, *minor),
                    _ => return Err(Error::InvalidTopology("a torus meeting off a quadric")),
                };
            let with = first.with == d0;
            let t0 = angle_of(if d0 { &first.pos[0] } else { &first.pos[1] });
            let t1 = angle_of(if dl { &last.pos[1] } else { &last.pos[0] });
            let sweep = if e.ends.is_none() {
                TAU
            } else {
                let s = if with { t1 - t0 } else { t0 - t1 };
                let s = s.rem_euclid(TAU);
                if s == 0.0 {
                    TAU
                } else {
                    s
                }
            };
            let sweep = if with { sweep } else { -sweep };
            let start = if e.ends.is_none() { 0.0 } else { t0 };
            if c.coaxial && !c.over_v {
                // A circle about the torus's axis at its root's height.
                let x = c
                    .at(&[int(1), crate::solid::split::zero()])
                    .ok_or(Error::ComputationLimit("a torus meeting's circle"))?;
                let cm = super::chain::curve_model(arr, first.curve, c.carrier, |_| true);
                let l = cm.f.local_q(&x);
                let (lu, lv, lw) = (l[0].to_f64(), l[1].to_f64(), l[2].to_f64());
                let circle = Frame3::new(
                    frame.point(Point2::default(), lw),
                    frame.normal(),
                    frame.x(),
                    cm.tolerance,
                )?;
                let radius = lu.hypot(lv);
                return Ok(if e.ends.is_none() {
                    Curve3::Circle {
                        frame: circle,
                        radius,
                    }
                } else {
                    Curve3::CircularArc {
                        frame: circle,
                        radius,
                        start_angle: start,
                        sweep_angle: sweep,
                    }
                });
            }
            Ok(Curve3::Toric(Box::new(crate::topology::Toric {
                frame: *frame,
                major: *major,
                minor: *minor,
                other,
                other_radius,
                other_sphere,
                other_half_angle,
                other_minor,
                over_v: c.over_v,
                window: c.window_angles(),
                start,
                sweep,
            })))
        }
        Crv::Torus(c) => {
            // S8d.3's spiric section on the stored torus (S9d.4a): the
            // plane in its frame, a unit normal; over its parameter's angle.
            // The torus's primitive model (a given face's view, S9e.3a).
            let m = super::chain::curve_model(arr, first.curve, c.carrier, |v| v.ring.is_some());
            let stored = &m
                .faces
                .iter()
                .find(|f| matches!(f.surf, Surf::Torus))
                .expect("a torus's patch")
                .stored;
            let Surface::Torus {
                frame,
                major,
                minor,
            } = stored
            else {
                return Err(Error::InvalidTopology("a torus section off a torus"));
            };
            let pl = c.plane.clone().map(|x| rational_f64(&x));
            let norm = (pl[0] * pl[0] + pl[1] * pl[1] + pl[2] * pl[2]).sqrt();
            let with = first.with == d0;
            let t0 = angle_of(if d0 { &first.pos[0] } else { &first.pos[1] });
            let t1 = angle_of(if dl { &last.pos[1] } else { &last.pos[0] });
            let sweep = if e.ends.is_none() {
                TAU
            } else {
                let s = if with { t1 - t0 } else { t0 - t1 };
                let s = s.rem_euclid(TAU);
                if s == 0.0 {
                    TAU
                } else {
                    s
                }
            };
            Ok(Curve3::Section(Box::new(crate::topology::Spiric {
                frame: *frame,
                major: *major,
                minor: *minor,
                plane: pl.map(|x| x / norm),
                over_v: c.over_v,
                sign: if c.plus { 1.0 } else { -1.0 },
                start: if e.ends.is_none() { 0.0 } else { t0 },
                sweep: if with { sweep } else { -sweep },
            })))
        }
        Crv::Circle(c) => {
            // A circle of a surd radius (S9d.1) on its basis's frame.
            let (Pos::Ang(p0), Pos::Ang(p1)) = (
                if d0 { &first.pos[0] } else { &first.pos[1] },
                if dl { &last.pos[1] } else { &last.pos[0] },
            ) else {
                unreachable!("a circle's places")
            };
            let with = first.with == d0;
            // A given result's circle (S9e.3a): on its stored circle.
            if let CurveRef::Edge(o, ei) = first.curve {
                if let Some(g) = &arr.models[o].given {
                    if let Some(c) = given_arc(g.curves[ei].as_ref(), e, points, with != g.flip[ei])
                    {
                        return Ok(c);
                    }
                }
            }
            let (t0, t1) = (c.angle(p0), c.angle(p1));
            let sweep = if e.ends.is_none() {
                TAU
            } else {
                let s = if with { t1 - t0 } else { t0 - t1 };
                let s = s.rem_euclid(TAU);
                if s == 0.0 {
                    TAU
                } else {
                    s
                }
            };
            let fl = |x: &V| {
                let a = x.clone().map(|y| rational_f64(&y));
                Vec3::new(a[0], a[1], a[2])
            };
            let (x, y) = (fl(&c.x), fl(&c.y));
            let o = fl(&c.c);
            let frame = Frame3::new(
                Point3::new(o.x, o.y, o.z),
                x.cross(y),
                x,
                arr.models[0].tolerance,
            )?;
            let radius = rational_f64(&c.r2).sqrt();
            Ok(if e.ends.is_none() {
                Curve3::Circle { frame, radius }
            } else {
                Curve3::CircularArc {
                    frame,
                    radius,
                    start_angle: t0,
                    sweep_angle: if with { sweep } else { -sweep },
                }
            })
        }
        Crv::Meet(m) => {
            // On the carrier's and the other's stored cylinders (S9c.2).
            let CurveRef::Section(si, _) = first.curve else {
                unreachable!("a meeting is a section")
            };
            let s = &arr.secs[si];
            let stored = [
                &arr.models[0].faces[s.fa].stored,
                &arr.models[1].faces[s.fb].stored,
            ];
            // A cylinder or cone (S9d.3b): its frame, radius and half angle.
            let ruled = |x: &Surface| match x {
                Surface::Cylinder { frame, radius } => Ok((*frame, *radius, 0.0)),
                Surface::Cone {
                    frame,
                    radius,
                    half_angle,
                } => Ok((*frame, *radius, *half_angle)),
                _ => Err(Error::InvalidTopology("a meeting off a ruled surface")),
            };
            let (frame, radius, half_angle) = ruled(stored[m.carrier])?;
            let (other, other_radius, other_sphere, other_half_angle) = match stored[1 - m.carrier]
            {
                Surface::Sphere { frame, radius } => (*frame, *radius, true, 0.0),
                s => {
                    let (f, r, a) = ruled(s)?;
                    (f, r, false, a)
                }
            };
            let with = first.with == d0;
            let t0 = angle_of(if d0 { &first.pos[0] } else { &first.pos[1] });
            let t1 = angle_of(if dl { &last.pos[1] } else { &last.pos[0] });
            let sweep = if e.ends.is_none() {
                TAU
            } else {
                let s = if with { t1 - t0 } else { t0 - t1 };
                let s = s.rem_euclid(TAU);
                if s == 0.0 {
                    TAU
                } else {
                    s
                }
            };
            // The stored frame's angle of the model's: its x axis turned.
            let cm = super::chain::curve_model(arr, first.curve, m.carrier, |_| true);
            let base = frame
                .x()
                .dot(cm.frame.y())
                .atan2(frame.x().dot(cm.frame.x()));
            Ok(Curve3::Meet(Box::new(crate::topology::Meet {
                frame,
                radius,
                half_angle,
                other,
                other_radius,
                other_sphere,
                other_half_angle,
                sign: if m.plus { 1.0 } else { -1.0 },
                start: if e.ends.is_none() { 0.0 } else { t0 - base },
                sweep: if with { sweep } else { -sweep },
            })))
        }
        Crv::Line { .. } => {
            let [s, t] = e
                .ends
                .ok_or(Error::InvalidTopology("a line without ends"))?;
            Ok(Curve3::LineSegment {
                start: points[&s],
                end: points[&t],
            })
        }
        Crv::Conic { c, a, b } => {
            // Start and end angles, and the turn along the chain.
            let with = first.with == d0;
            let t0 = angle_of(if d0 { &first.pos[0] } else { &first.pos[1] });
            let t1 = angle_of(if dl { &last.pos[1] } else { &last.pos[0] });
            let sweep = if e.ends.is_none() {
                TAU
            } else {
                let s = if with { t1 - t0 } else { t0 - t1 };
                let s = s.rem_euclid(TAU);
                if s == 0.0 {
                    TAU
                } else {
                    s
                }
            };
            let sweep = if with { sweep } else { -sweep };
            let fl = |x: &V| x.clone().map(|y| rational_f64(&y));
            // A given result's edge (S9e.1): on its stored circle or ellipse
            // (turning its frame's way, S9e.2).
            if let CurveRef::Edge(o, ei) = first.curve {
                if let Some(g) = &arr.models[o].given {
                    if let Some(c) = given_arc(g.curves[ei].as_ref(), e, points, with != g.flip[ei])
                    {
                        return Ok(c);
                    }
                }
            }
            // A model arc: on its cap's arc frame.
            if let (CurveRef::Edge(o, ei), None) = (first.curve, &arr.models[first_op(first)].given)
            {
                let m = &arr.models[o];
                if let (EdgeKind::Rim(high, _), Some(fun)) = (m.edges[ei].kind, &m.funnel) {
                    // A cone's rim (S9d.3a): about its axis at its end.
                    let (r, h) = if high {
                        (&fun.t, rational_f64(&m.hi))
                    } else {
                        (&fun.b, 0.0)
                    };
                    let frame = Frame3::new(
                        m.frame.point(Point2::default(), h),
                        m.frame.normal(),
                        m.frame.x(),
                        m.tolerance,
                    )?;
                    let radius = rational_f64(r);
                    return Ok(if e.ends.is_none() {
                        Curve3::Circle { frame, radius }
                    } else {
                        Curve3::CircularArc {
                            frame,
                            radius,
                            start_angle: t0,
                            sweep_angle: sweep,
                        }
                    });
                }
                let EdgeKind::Cap(high, b2, j) = m.edges[ei].kind else {
                    unreachable!("an arc edge is a cap edge")
                };
                let Seg::Arc { c: ac, r, .. } = &m.bounds[b2].segs[j] else {
                    unreachable!("an arc segment")
                };
                let h = rational_f64(if high { &m.hi } else { &m.lo });
                let frame = Frame3::new(
                    m.frame
                        .point(Point2::new(rational_f64(&ac[0]), rational_f64(&ac[1])), h),
                    m.frame.normal(),
                    m.frame.x(),
                    m.tolerance,
                )?;
                let radius = rational_f64(r);
                return Ok(if e.ends.is_none() {
                    Curve3::Circle { frame, radius }
                } else {
                    Curve3::CircularArc {
                        frame,
                        radius,
                        start_angle: t0,
                        sweep_angle: sweep,
                    }
                });
            }
            // A section: principal axes of c + a cos t + b sin t.
            let (cf, af, bf) = (fl(c), fl(a), fl(b));
            let d = |x: &[f64; 3], y: &[f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
            let (aa, ab, bb) = (d(&af, &af), d(&af, &bf), d(&bf, &bf));
            let phi = 0.5 * (2.0 * ab).atan2(aa - bb);
            let (cp, sp) = (phi.cos(), phi.sin());
            let mut major = [0, 1, 2].map(|k| af[k] * cp + bf[k] * sp);
            let mut minor = [0, 1, 2].map(|k| -af[k] * sp + bf[k] * cp);
            let mut shift = phi;
            if d(&minor, &minor) > d(&major, &major) {
                let m2 = minor;
                minor = major.map(|x| -x);
                major = m2;
                shift += FRAC_PI_2;
            }
            let (lm, ln) = (d(&major, &major).sqrt(), d(&minor, &minor).sqrt());
            let x = Vec3::new(major[0], major[1], major[2]);
            let y = Vec3::new(minor[0], minor[1], minor[2]);
            let frame = Frame3::new(
                Point3::new(cf[0], cf[1], cf[2]),
                x.cross(y),
                x,
                arr.models[0].tolerance,
            )?;
            let _ = points;
            Ok(Curve3::EllipseArc {
                frame,
                major: lm,
                minor: ln,
                start_angle: if e.ends.is_none() { 0.0 } else { t0 - shift },
                sweep_angle: sweep,
            })
        }
    }
}

/// The operand of an arrangement edge's model edge (0 for a section).
fn first_op(g: &GEdge) -> usize {
    match g.curve {
        CurveRef::Edge(o, _) => o,
        CurveRef::Section(..) => 0,
    }
}

/// A given result's edge's piece on its stored circle or ellipse (the
/// frame the result rounded), its ends' angles measured on it from their
/// rounded points: none for another curve.
fn given_arc(
    curve: Option<&Curve3>,
    e: &REdge,
    points: &BTreeMap<usize, Point3>,
    with: bool,
) -> Option<Curve3> {
    let turn = |t0: f64, t1: f64| {
        let s = if with { t1 - t0 } else { t0 - t1 }.rem_euclid(TAU);
        let s = if s == 0.0 { TAU } else { s };
        if with {
            s
        } else {
            -s
        }
    };
    match curve? {
        Curve3::Circle { frame, radius } | Curve3::CircularArc { frame, radius, .. } => {
            let Some([s, t]) = e.ends else {
                return Some(Curve3::Circle {
                    frame: *frame,
                    radius: *radius,
                });
            };
            let angle = |p: Point3| {
                let [x, y, _] = frame.coordinates(p);
                y.atan2(x)
            };
            let (t0, t1) = (angle(points[&s]), angle(points[&t]));
            Some(Curve3::CircularArc {
                frame: *frame,
                radius: *radius,
                start_angle: t0,
                sweep_angle: turn(t0, t1),
            })
        }
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            ..
        } => {
            let Some([s, t]) = e.ends else {
                return Some(Curve3::EllipseArc {
                    frame: *frame,
                    major: *major,
                    minor: *minor,
                    start_angle: 0.0,
                    sweep_angle: if with { TAU } else { -TAU },
                });
            };
            let angle = |p: Point3| {
                let [x, y, _] = frame.coordinates(p);
                (y / minor).atan2(x / major)
            };
            let (t0, t1) = (angle(points[&s]), angle(points[&t]));
            Some(Curve3::EllipseArc {
                frame: *frame,
                major: *major,
                minor: *minor,
                start_angle: t0,
                sweep_angle: turn(t0, t1),
            })
        }
        _ => None,
    }
}

/// A loop's fins: each use's pcurve on the face's stored surface, lifted
/// continuously along the loop on a cylinder.
fn loop_fins(
    p: &TopologyParts,
    uses: &[(usize, bool)],
    edge_id: &BTreeMap<usize, EdgeId>,
    surface: &Surface,
) -> Result<Vec<Fin>> {
    let mut out = Vec::new();
    let mut lift: Option<Point2> = None;
    for &(ri, fwd) in uses {
        let edge = edge_id[&ri];
        let curve = &p.edges[edge.0].curve;
        let sense = if fwd {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        let pcurve = match (surface, curve) {
            (Surface::Plane(frame), Curve3::EllipseArc { .. }) => {
                let reversed = !fwd;
                let start = curve.point(if reversed { 1.0 } else { 0.0 });
                let [x, y, _] = frame.coordinates(start);
                Curve2::Projection(Box::new(
                    Projection::new(
                        curve.clone(),
                        surface.clone(),
                        reversed,
                        Point2::new(x, y),
                        16,
                    )
                    .ok_or(Error::PrecisionLoss)?,
                ))
            }
            (Surface::Plane(frame), _) => plane_pcurve(curve, sense, *frame),
            _ => {
                let reversed = !fwd;
                let start = curve.point(if reversed { 1.0 } else { 0.0 });
                let mut uv = Projection::inverse(surface, start).ok_or(Error::PrecisionLoss)?;
                if let Some(prev) = lift {
                    uv.x += TAU * ((prev.x - uv.x) / TAU).round();
                    // A torus's v is periodic too (S9d.4a).
                    if matches!(surface, Surface::Torus { .. }) {
                        uv.y += TAU * ((prev.y - uv.y) / TAU).round();
                    }
                }
                let pc = match cylinder_pcurve(surface, curve, reversed, uv)
                    .or_else(|| sphere_pcurve(surface, curve, reversed, lift.map(|l| l.x)))
                    .or_else(|| cone_pcurve(surface, curve, reversed, uv))
                {
                    Some(pc) => pc,
                    // More anchors where 16 do not pin the lift (a curve
                    // passing near a sphere's pole, S9d.4b.2).
                    None => Curve2::Projection(Box::new(
                        [16, 64, 256]
                            .into_iter()
                            .find_map(|n| {
                                Projection::new(curve.clone(), surface.clone(), reversed, uv, n)
                            })
                            .ok_or(Error::PrecisionLoss)?,
                    )),
                };
                lift = Some(pc.point(1.0));
                pc
            }
        };
        out.push(Fin {
            edge,
            sense,
            pcurve,
            enclosure: None,
        });
    }
    Ok(out)
}

/// A cylinder face's pcurve of an edge from `uv` (its start's lifted
/// parameters, along the use): a line for a generatrix or a circle about
/// the axis, the sinusoid `v = a0 + a1 cos u + a2 sin u` for a plane's
/// section whose parameter is the cylinder's angle (checked along it);
/// `None` otherwise (an exact projection instead).
fn cylinder_pcurve(
    surface: &Surface,
    curve: &Curve3,
    reversed: bool,
    uv: Point2,
) -> Option<Curve2> {
    let Surface::Cylinder { frame, radius } = surface else {
        return None;
    };
    let at = |f: f64| curve.point(if reversed { 1.0 - f } else { f });
    let end = Projection::inverse(surface, at(1.0))?;
    let (plane, sweep) = match curve {
        Curve3::LineSegment { .. } => {
            // A generatrix: u constant.
            return ((end.x - uv.x)
                .rem_euclid(TAU)
                .min(TAU - (end.x - uv.x).rem_euclid(TAU))
                < 1e-9)
                .then_some(Curve2::LineSegment {
                    start: uv,
                    end: Point2::new(uv.x, end.y),
                });
        }
        Curve3::Circle { frame: f, .. } => (*f, TAU),
        Curve3::CircularArc {
            frame: f,
            sweep_angle,
            ..
        }
        | Curve3::EllipseArc {
            frame: f,
            sweep_angle,
            ..
        } => (*f, *sweep_angle),
        _ => return None,
    };
    let sweep = if reversed { -sweep } else { sweep };
    // The curve's plane m . (P - p) = 0 on the cylinder's points.
    let (m, p0) = (plane.normal(), plane.origin());
    let (o, x, y, n) = (frame.origin(), frame.x(), frame.y(), frame.normal());
    let mn = m.dot(n);
    if mn.abs() < 1e-9 {
        return None;
    }
    let a = [
        m.dot(p0 - o) / mn,
        -radius * m.dot(x) / mn,
        -radius * m.dot(y) / mn,
    ];
    let pc = if a[1].abs() <= 1e-15 * a[0].abs().max(*radius)
        && a[2].abs() <= 1e-15 * a[0].abs().max(*radius)
    {
        Curve2::LineSegment {
            start: Point2::new(uv.x, a[0]),
            end: Point2::new(uv.x + sweep, a[0]),
        }
    } else {
        Curve2::Sinusoid {
            start: uv.x,
            sweep,
            a,
        }
    };
    // The same points at the same fractions.
    let tol = 1e-9 * (1.0 + radius + o.to_array().iter().fold(0.0f64, |m, v| m.max(v.abs())));
    for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let q = pc.point(f);
        let s = frame.point(Point2::new(radius * q.x.cos(), radius * q.x.sin()), q.y);
        if (s - at(f)).length() > tol {
            return None;
        }
    }
    Some(pc)
}

/// A sphere face's pcurve of an edge (S9d.1): a meridian's `u` constant
/// (taken inside it, so an end at a pole needs none), a parallel's `v`
/// constant, lines checked against the edge at the same fractions; `None`
/// otherwise (an exact projection instead). `u` is lifted near `prev`.
fn sphere_pcurve(
    surface: &Surface,
    curve: &Curve3,
    reversed: bool,
    prev: Option<f64>,
) -> Option<Curve2> {
    let Surface::Sphere { frame, radius } = surface else {
        return None;
    };
    if !matches!(curve, Curve3::Circle { .. } | Curve3::CircularArc { .. }) {
        return None;
    }
    let at = |f: f64| curve.point(if reversed { 1.0 - f } else { f });
    let fr = [0.0, 0.25, 0.5, 0.75, 1.0];
    let uv: Vec<Point2> = fr
        .iter()
        .map(|&f| Projection::inverse(surface, at(f)))
        .collect::<Option<_>>()?;
    let near = |x: f64, target: f64| x + TAU * ((target - x) / TAU).round();
    let tol = 1e-9
        * (1.0
            + radius
            + frame
                .origin()
                .to_array()
                .iter()
                .fold(0.0f64, |m, v| m.max(v.abs())));
    let check = |pc: &Curve2| {
        fr.iter().all(|&f| {
            let q = pc.point(f);
            let s = frame.point(
                Point2::new(
                    radius * q.y.cos() * q.x.cos(),
                    radius * q.y.cos() * q.x.sin(),
                ),
                radius * q.y.sin(),
            );
            (s - at(f)).length() <= tol
        })
    };
    // A meridian: u equal at the interior samples.
    let u_in = [uv[1].x, uv[2].x, uv[3].x];
    let same_u = u_in
        .iter()
        .all(|u| (near(*u, u_in[0]) - u_in[0]).abs() < 1e-9);
    if same_u {
        let u = prev.map_or(u_in[0], |p| near(u_in[0], p));
        let pc = Curve2::LineSegment {
            start: Point2::new(u, uv[0].y),
            end: Point2::new(u, uv[4].y),
        };
        if check(&pc) {
            return Some(pc);
        }
    }
    // A parallel: v equal throughout, u turning with the curve.
    if uv.iter().all(|p| (p.y - uv[0].y).abs() < 1e-9) && uv[0].y.cos() > 1e-6 {
        let mut u = prev.map_or(uv[0].x, |p| near(uv[0].x, p));
        let u0 = u;
        for w in uv.windows(2) {
            let d = near(w[1].x, w[0].x) - w[0].x;
            u += d;
        }
        let pc = Curve2::LineSegment {
            start: Point2::new(u0, uv[0].y),
            end: Point2::new(u, uv[0].y),
        };
        if check(&pc) {
            return Some(pc);
        }
    }
    None
}

/// A cone face's pcurve of a circle about its axis (S9d.3a): `v`
/// constant, `u` turning with the curve from `uv` (its start's lifted
/// parameters), checked at the same fractions; `None` otherwise.
fn cone_pcurve(surface: &Surface, curve: &Curve3, reversed: bool, uv: Point2) -> Option<Curve2> {
    let Surface::Cone { frame, .. } = surface else {
        return None;
    };
    if !matches!(curve, Curve3::Circle { .. } | Curve3::CircularArc { .. }) {
        return None;
    }
    let at = |f: f64| curve.point(if reversed { 1.0 - f } else { f });
    let fr = [0.0, 0.25, 0.5, 0.75, 1.0];
    let uvs: Vec<Point2> = fr
        .iter()
        .map(|&f| Projection::inverse(surface, at(f)))
        .collect::<Option<_>>()?;
    let scale = 1.0
        + frame
            .origin()
            .to_array()
            .iter()
            .fold(0.0f64, |m, v| m.max(v.abs()));
    if uvs.iter().any(|p| (p.y - uvs[0].y).abs() > 1e-9 * scale) {
        return None;
    }
    let near = |x: f64, target: f64| x + TAU * ((target - x) / TAU).round();
    let mut u = uv.x;
    for w in uvs.windows(2) {
        u += near(w[1].x, w[0].x) - w[0].x;
    }
    let pc = Curve2::LineSegment {
        start: Point2::new(uv.x, uvs[0].y),
        end: Point2::new(u, uvs[0].y),
    };
    fr.iter()
        .all(|&f| (surface.point(pc.point(f)) - at(f)).length() <= 1e-9 * scale)
        .then_some(pc)
}

/// A plane's section of a cone rounded once (S9d.3a): S8d.2's ellipse,
/// hyperbola or parabola (a circle about the axis for a plane normal to
/// it), its parameter running as the chain runs.
fn cone_curve3(
    arr: &Arr,
    m: &super::model::Prism,
    e: &REdge,
    c: &super::cone::ConeSec,
    points: &BTreeMap<usize, Point3>,
) -> Result<Curve3> {
    let fun = m.funnel.as_ref().expect("a cone's section on a cone");
    let (g0, d0) = e.parts[0];
    // The chain's direction at its start, from its first part's samples.
    let pts = arr.samples((g0, d0));
    let (p0, p1) = (
        Point3::new(pts[0][0], pts[0][1], pts[0][2]),
        Point3::new(pts[1][0], pts[1][1], pts[1][2]),
    );
    let ends = e.ends.map(|[s, t]| (points[&s], points[&t]));
    if c.normal_to_axis() {
        // A circle about the axis.
        let (rho, w) = c.circle();
        let frame = Frame3::new(
            m.frame.point(Point2::default(), rational_f64(&w)),
            m.frame.normal(),
            m.frame.x(),
            m.tolerance,
        )?;
        let radius = rational_f64(&rho);
        let angle = |p: Point3| {
            let [x, y, _] = frame.coordinates(p);
            y.atan2(x)
        };
        let turn =
            (angle(p1) - angle(p0) + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI;
        let dir = if turn >= 0.0 { 1.0 } else { -1.0 };
        return Ok(match ends {
            None => Curve3::Circle { frame, radius },
            Some((s, t)) => {
                let (t0, t1) = (angle(s), angle(t));
                let mut sweep = if dir > 0.0 { t1 - t0 } else { t0 - t1 }.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                Curve3::CircularArc {
                    frame,
                    radius,
                    start_angle: t0,
                    sweep_angle: dir * sweep,
                }
            }
        });
    }
    let (section, _) = crate::solid::split::conic::cone_conic(
        m.frame,
        m.tolerance,
        rational_f64(&fun.b),
        rational_f64(&fun.t),
        rational_f64(&m.hi),
        &c.plane,
    )?;
    if matches!(section, crate::solid::split::conic::Section::Rulings) {
        return Err(Error::Degenerate(
            "a plane within the resolution of a cone's apex",
        ));
    }
    let periodic = section.periodic();
    let wrap = |d: f64| {
        if periodic {
            (d + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI
        } else {
            d
        }
    };
    let dir = if wrap(section.param(p1) - section.param(p0)) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Ok(match ends {
        None => section.arc(section.param(p0), dir * TAU),
        Some((s, t)) => {
            let (t0, t1) = (section.param(s), section.param(t));
            let sweep = if periodic {
                let mut x = if dir > 0.0 { t1 - t0 } else { t0 - t1 }.rem_euclid(TAU);
                if x == 0.0 {
                    x = TAU;
                }
                dir * x
            } else {
                t1 - t0
            };
            section.arc(t0, sweep)
        }
    })
}

/// A torus face's loops on one sheet of the cover (S9d.4a): with a loop
/// winding in `u` as reference, each other loop's `v` shifted by whole
/// turns into the turn above it (a `+u` reference on a forward face: its
/// material above) or below it; with only loops winding in `v`, likewise
/// in `u` (a `-v` reference: material towards `+u`).
fn torus_sheets(p: &mut TopologyParts, loops: &[LoopId], reversed: bool) {
    let info = |p: &TopologyParts, l: LoopId| -> Option<([i32; 2], [f64; 2])> {
        let Loop::Edges { fins, winding } = &p.loops[l.0] else {
            return None;
        };
        let pts: Vec<Point2> = fins
            .iter()
            .flat_map(|f| (0..8).map(move |i| (f, i as f64 / 8.0)))
            .map(|(f, t)| p.fins[f.0].pcurve.point(t))
            .collect();
        let n = pts.len() as f64;
        Some((
            *winding,
            [
                pts.iter().map(|q| q.x).sum::<f64>() / n,
                pts.iter().map(|q| q.y).sum::<f64>() / n,
            ],
        ))
    };
    let all: Vec<(LoopId, [i32; 2], [f64; 2])> = loops
        .iter()
        .filter_map(|&l| info(p, l).map(|(w, m)| (l, w, m)))
        .collect();
    // (coordinate, reference loop, material above it)
    let reference = all
        .iter()
        .find(|x| x.1[0] != 0)
        .map(|x| (1usize, x.0, x.1[0] > 0, x.2[1]))
        .or_else(|| {
            all.iter()
                .find(|x| x.1[1] != 0)
                .map(|x| (0usize, x.0, x.1[1] < 0, x.2[0]))
        });
    let Some((axis, rid, above, rmean)) = reference else {
        return;
    };
    for (l, _, mean) in &all {
        if *l == rid {
            continue;
        }
        let m = mean[axis];
        // A reversed face's material lies on the loops' other side.
        let target = if above != reversed {
            rmean + std::f64::consts::PI
        } else {
            rmean - std::f64::consts::PI
        };
        let turns = ((target - m) / TAU).round();
        if turns == 0.0 {
            continue;
        }
        let Loop::Edges { fins, .. } = &p.loops[l.0] else {
            continue;
        };
        let k = TAU * turns;
        for f in fins.clone() {
            let shift = |q: &mut Point2| {
                if axis == 0 {
                    q.x += k;
                } else {
                    q.y += k;
                }
            };
            match &mut p.fins[f.0].pcurve {
                Curve2::Projection(pr) => pr.lifts.iter_mut().for_each(shift),
                Curve2::LineSegment { start, end } => {
                    shift(start);
                    shift(end);
                }
                _ => {}
            }
        }
    }
}

/// A torus face's loops none of which winds (S9d.4b.1: a wedge's face over
/// more than half a turn, its holes' pcurves lifted from the surface's
/// principal angles): each other loop shifted by whole turns in `u` and `v`
/// to lie inside the first, its outer loop, on the cover.
fn torus_holes_in(p: &mut TopologyParts, loops: &[LoopId]) {
    let samples = |p: &TopologyParts, l: LoopId| -> Option<Vec<[f64; 2]>> {
        let Loop::Edges { fins, winding } = &p.loops[l.0] else {
            return None;
        };
        (*winding == [0, 0]).then(|| {
            fins.iter()
                .flat_map(|f| (0..8).map(move |i| (f, i as f64 / 8.0)))
                .map(|(f, t)| {
                    let q = p.fins[f.0].pcurve.point(t);
                    [q.x, q.y]
                })
                .collect()
        })
    };
    let all: Option<Vec<Vec<[f64; 2]>>> = loops.iter().map(|&l| samples(p, l)).collect();
    let Some(all) = all else {
        return;
    };
    let Some(outer) = all.first() else {
        return;
    };
    for (k, &l) in loops.iter().enumerate().skip(1) {
        let probe = all[k][0];
        let shift = (-2..=2i32)
            .flat_map(|a| (-2..=2i32).map(move |b| (a, b)))
            .find(|&(a, b)| {
                let x = [probe[0] + TAU * f64::from(a), probe[1] + TAU * f64::from(b)];
                super::graph::winding(outer, x).0
            });
        let Some((a, b)) = shift.filter(|s| *s != (0, 0)) else {
            continue;
        };
        let (du, dv) = (TAU * f64::from(a), TAU * f64::from(b));
        let Loop::Edges { fins, .. } = &p.loops[l.0] else {
            continue;
        };
        for f in fins.clone() {
            let shift = |q: &mut Point2| {
                q.x += du;
                q.y += dv;
            };
            match &mut p.fins[f.0].pcurve {
                Curve2::Projection(pr) => pr.lifts.iter_mut().for_each(shift),
                Curve2::LineSegment { start, end } => {
                    shift(start);
                    shift(end);
                }
                _ => {}
            }
        }
    }
}

/// A torus loop's turns about the tube: its pcurves' change of `v`.
fn turns_v(p: &TopologyParts, fins: &[FinId]) -> i32 {
    let mut total = 0.0;
    for f in fins {
        let pc = &p.fins[f.0].pcurve;
        total += pc.point(1.0).y - pc.point(0.0).y;
    }
    (total / TAU).round() as i32
}

/// A cylinder loop's turns about the axis: its pcurves' change of `u`,
/// with the jumps between them (a sphere's or cone's pole, where `u` turns
/// freely: S9d.4c, a loop through a pole whose pcurves alone turn half a
/// turn), the last pcurve's end against the first's start.
fn turns(p: &TopologyParts, fins: &[FinId], _surface: &Surface) -> i32 {
    let (Some(first), Some(last)) = (fins.first(), fins.last()) else {
        return 0;
    };
    let total = p.fins[last.0].pcurve.point(1.0).x - p.fins[first.0].pcurve.point(0.0).x;
    (total / TAU).round() as i32
}

impl Arr {
    /// A face's loops from half-edges: at a vertex with several ways on,
    /// the first clockwise from the way back.
    pub(super) fn relink(
        &self,
        o: usize,
        hs: &[(usize, bool)],
        face_of: &BTreeMap<(usize, bool), usize>,
    ) -> Result<Vec<Vec<(usize, bool)>>> {
        let start = |h: (usize, bool)| self.edges[h.0].ends[if h.1 { 0 } else { 1 }];
        let end = |h: (usize, bool)| self.edges[h.0].ends[if h.1 { 1 } else { 0 }];
        let mut out_of: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
        for &h in hs {
            out_of.entry(start(h)).or_default().push(h);
        }
        let mut used: BTreeSet<(usize, bool)> = BTreeSet::new();
        let mut loops = Vec::new();
        for &h0 in hs {
            if used.contains(&h0) {
                continue;
            }
            used.insert(h0);
            let mut lp = vec![h0];
            let mut h = h0;
            loop {
                let v = end(h);
                let outs = out_of.get(&v).map_or(&[][..], |x| x);
                let next = if outs.len() == 1 {
                    outs[0]
                } else {
                    self.next_on(o, face_of[&h], v, h, outs)?
                };
                if next == h0 {
                    break;
                }
                if !used.insert(next) {
                    return Err(Error::InvalidTopology("a joined face's loops do not close"));
                }
                lp.push(next);
                h = next;
            }
            loops.push(lp);
        }
        let _ = Ordering::Equal;
        Ok(loops)
    }
}
