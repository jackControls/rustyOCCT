//! S9a.2: a Boolean of two prisms in one frame whose result is a stack of
//! slabs of different regions, built as general bodies.
//!
//! The heights of both prisms' ends cut the result into slabs; in each the
//! result's region is a set function of the two profiles (`A ∪ B`, `A`,
//! `B`, `A - B` or nothing), so one arrangement of both profiles' boundaries
//! serves every slab. Each arrangement piece knows whether each side lies
//! in `A` and in `B`; a piece is a wall in a slab where the result holds one
//! side and not the other, and the result's horizontal faces at a height are
//! the regions it holds on one side of the height only, traced from the
//! same arrangement. Walls on one line or circle facing one way join across
//! slabs and pieces where nothing else meets them; an edge keeps no vertex
//! where it runs straight on (collinear lines, arcs of one circle, a
//! vertical line) between the same two faces: the unified result, as the
//! kernel's single prisms are.
//!
//! Names follow provenance as in S9a.1 (`solid/boolean.rs`): each face,
//! edge and vertex continues the input entities it is a part of (walls on
//! an input's wall facing its way, caps on an input's cap, edges and
//! vertices on the inputs' cap edges, vertical edges and vertices) and lies
//! on the others (a cut's tool, which faces the other way, and an input's
//! wall or cap an entity runs inside).
use super::{At, What};
use crate::identity::OperationId;
use crate::identity::{EntityKind, Role};
use crate::profile::boolean::{arrange, boolean, select_with, Curve, Op2, Operand, PId, PieceView};
use crate::solid::split::q;
use crate::solid::{Construction, MassProperties, Solid};
use crate::topology::Topology;
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side, Slot, Surface, TopologyParts,
    Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Profile, Result, Tolerance};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

/// An input entity: the operand's entity at (end, what, boundary, element).
pub(super) type Key = (Operand, At, What, usize, usize);

/// A result slot's plan: the input entities it continues and those it lies
/// on, its kind and role.
pub(super) type KeyPlan = (Slot, Vec<Key>, Vec<Key>, EntityKind, Role);

/// One connected solid of a stack: its parts, each slot's provenance and
/// its height range.
pub(super) struct Component {
    pub(super) parts: TopologyParts,
    pub(super) plans: Vec<KeyPlan>,
    pub(super) heights: [f64; 2],
}

/// A stack's inputs (the tool's profile translated into the object's
/// frame), enough to build it again in another frame; `index` is the
/// component.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Stack {
    pub(super) a: Profile,
    pub(super) b: Profile,
    pub(super) ha: [f64; 2],
    pub(super) hb: [f64; 2],
    pub(super) op: Op2,
    pub(super) index: usize,
}

impl Stack {
    pub(crate) fn tolerance(&self) -> Tolerance {
        self.a.tolerance()
    }

    /// The same component of the same stack in another frame, its ids
    /// external (the caller restores them), with its mass properties when
    /// known (a rigid motion's: its source's, moved).
    pub(crate) fn rebuilt_with(
        &self,
        operation: OperationId,
        frame: Frame3,
        known: Option<MassProperties>,
    ) -> Result<Solid> {
        let mut components = build(self, frame)?;
        if self.index >= components.len() {
            return Err(Error::InvalidTopology("a stack rebuilt differently"));
        }
        let component = components.swap_remove(self.index);
        self.solid(component, frame, operation, known)
    }

    /// A component as a solid (external ids).
    pub(super) fn solid(
        &self,
        component: Component,
        frame: Frame3,
        operation: OperationId,
        known: Option<MassProperties>,
    ) -> Result<Solid> {
        let topology =
            Topology::from_parts(component.parts.with_measured_enclosures(), self.tolerance())
                .map_err(|issues| {
                    Error::InvalidTopology(issues.first().map_or("a stack", |i| i.kind.name()))
                })?;
        let mass = match known {
            Some(mass) => mass,
            None => topology
                .mass_enclosure()
                .ok_or(Error::Unrepresentable("a stack's mass properties"))?
                .midpoints(),
        };
        let bounds = crate::solid::split::edge_bounds(&topology);
        Ok(Solid {
            construction: Construction::Stack(Box::new(self.clone())),
            frame,
            start: component.heights[0],
            end: component.heights[1],
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// The result's material at a slab: whether a point in `A` (`m[0]`) and
    /// in `B` (`m[1]`) lies in it; nothing outside the slabs.
    fn holds(&self, heights: &[f64], s: isize, m: [bool; 2]) -> bool {
        if s < 0 || s as usize + 1 >= heights.len() {
            return false;
        }
        let s = s as usize;
        let present = |h: [f64; 2]| h[0] <= heights[s] && heights[s + 1] <= h[1];
        let (a, b) = (m[0] && present(self.ha), m[1] && present(self.hb));
        match self.op {
            Op2::Fuse => a || b,
            Op2::Cut => a && !b,
            Op2::Common => a && b,
        }
    }

    fn heights(&self) -> Vec<f64> {
        let mut h = vec![self.ha[0], self.ha[1], self.hb[0], self.hb[1]];
        h.sort_by(f64::total_cmp);
        h.dedup();
        h
    }

    /// Inside where the result's slab at the point's height holds it; on
    /// the boundary within the resolution of a face (`range`: the
    /// component's heights).
    pub(crate) fn classify(
        &self,
        [u, v, w]: [f64; 3],
        range: [f64; 2],
        tolerance: Tolerance,
    ) -> Result<crate::Location> {
        use crate::decide::sum_le;
        use crate::Location;
        for x in [u, v, w] {
            crate::math::finite(x, "coordinate")?;
        }
        let tol = tolerance.linear();
        if !(sum_le(&[range[0], -tol], &[w]) && sum_le(&[w], &[range[1], tol])) {
            return Ok(Location::Outside);
        }
        let heights = self.heights();
        let p = Point2::new(u, v);
        // Each profile's possible memberships: both on its boundary.
        let options = |l: Location| match l {
            Location::Inside => vec![true],
            Location::Outside => vec![false],
            Location::Boundary => vec![false, true],
        };
        let (la, lb) = (options(self.a.classify(p)?), options(self.b.classify(p)?));
        let mut seen = BTreeSet::new();
        // The component's slabs within the resolution of the height, and
        // nothing beyond its ends.
        for s in 0..heights.len() - 1 {
            let (lo, hi) = (heights[s].max(range[0]), heights[s + 1].min(range[1]));
            if lo >= hi || !(sum_le(&[lo, -tol], &[w]) && sum_le(&[w], &[hi, tol])) {
                continue;
            }
            for a in &la {
                for b in &lb {
                    seen.insert(self.holds(&heights, s as isize, [*a, *b]));
                }
            }
        }
        let near = |h: f64| sum_le(&[w, -h], &[tol]) && sum_le(&[h, -w], &[tol]);
        if near(range[0]) || near(range[1]) {
            seen.insert(false);
        }
        Ok(match (seen.contains(&true), seen.contains(&false)) {
            (true, false) => Location::Inside,
            (true, true) => Location::Boundary,
            _ => Location::Outside,
        })
    }
}

/// A fine edge: a piece at a level, or a point's vertical over a slab.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Fine {
    H(usize, usize),
    V(PId, usize),
}

/// A face cell: a piece's wall over a slab, or a horizontal face at a
/// level (upward, its index among that level's).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Cell {
    Wall(usize, usize),
    Cap(usize, bool, usize),
}

/// A point at a level.
type Node = (PId, usize);

/// A horizontal face's region and its region-left cycles of pieces, each
/// with whether it runs along its segment.
type CapFace = (Profile, Vec<Vec<(usize, bool)>>);

fn find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut y = x;
    while parent[y] != r {
        let next = parent[y];
        parent[y] = r;
        y = next;
    }
    r
}

fn union(parent: &mut [usize], x: usize, y: usize) {
    let (a, b) = (find(parent, x), find(parent, y));
    if a != b {
        parent[a.max(b)] = a.min(b);
    }
}

/// A piece's circle, if on one.
fn circle(v: &PieceView) -> Option<(Point2, f64)> {
    match v.curve {
        Curve::Arc { center, radius, .. } | Curve::Circle { center, radius } => {
            Some((center, radius))
        }
        Curve::Line => None,
    }
}

/// Whether a piece runs counter-clockwise along its stored direction.
fn stored_ccw(v: &PieceView) -> bool {
    match v.curve {
        Curve::Arc { ccw, .. } => ccw,
        _ => true,
    }
}

/// Whether two pieces lie on one line or one circle.
fn same_carrier(x: &PieceView, y: &PieceView) -> bool {
    match (circle(x), circle(y)) {
        (None, None) => {
            let d = |v: &PieceView| [q(v.e.x) - q(v.p.x), q(v.e.y) - q(v.p.y)];
            let (a, b) = (d(x), d(y));
            let cross = &a[0] * &b[1] - &a[1] * &b[0];
            // Collinear directions through a shared point: one line.
            cross == q(0.0) && {
                let w = [q(y.p.x) - q(x.p.x), q(y.p.y) - q(x.p.y)];
                &a[0] * &w[1] - &a[1] * &w[0] == q(0.0)
            }
        }
        (Some(c), Some(k)) => c == k,
        _ => false,
    }
}

/// Builds a stack's components in `frame`.
#[allow(clippy::too_many_lines)]
pub(super) fn build(stack: &Stack, frame: Frame3) -> Result<Vec<Component>> {
    let tolerance = stack.tolerance();
    let tol = tolerance.linear();
    let heights = stack.heights();
    let slabs = heights.len() - 1;
    for s in 0..slabs {
        if heights[s + 1] - heights[s] <= tol {
            return Err(Error::Degenerate("a slab thinner than the resolution"));
        }
    }
    let arranged = arrange(&[&stack.a], &[&stack.b], tolerance)?;
    let views = arranged.pieces();
    let holds = |s: isize, m: [bool; 2]| stack.holds(&heights, s, m);
    let n = frame.normal();
    let at = |p: PId, k: usize| frame.point(arranged.position(p), heights[k]);

    // Cells: walls, then each level's horizontal faces.
    let mut cells: Vec<Cell> = Vec::new();
    let mut material_left: BTreeMap<(usize, usize), bool> = BTreeMap::new();
    for (g, v) in views.iter().enumerate() {
        if !v.representative {
            continue;
        }
        for s in 0..slabs {
            let (l, r) = (holds(s as isize, v.left), holds(s as isize, v.right));
            if l != r {
                cells.push(Cell::Wall(g, s));
                material_left.insert((g, s), l);
            }
        }
    }
    // Per level and way (up: material below only), its faces' profiles
    // and region-left cycles.
    let mut caps: BTreeMap<(usize, bool), Vec<CapFace>> = BTreeMap::new();
    for k in 0..=slabs {
        for up in [true, false] {
            let way = |a: bool, b: bool| {
                let (below, above) = (holds(k as isize - 1, [a, b]), holds(k as isize, [a, b]));
                if up {
                    below && !above
                } else {
                    above && !below
                }
            };
            let any = [(false, false), (true, false), (false, true), (true, true)]
                .iter()
                .any(|&(a, b)| way(a, b));
            if !any {
                continue;
            }
            let traced = select_with(&arranged, &way, false, tolerance)?;
            let faces: Vec<_> = traced
                .pieces
                .into_iter()
                .map(|p| (p.profile, p.cycles))
                .collect();
            for i in 0..faces.len() {
                cells.push(Cell::Cap(k, up, i));
            }
            caps.insert((k, up), faces);
        }
    }

    // Each cell's half-edges: fine edges with whether they run along them
    // (a piece's stored direction, upward), the cell on their left seen
    // from outside the material.
    let ends = |g: usize| views[g].from.zip(views[g].to);
    let mut half: Vec<Vec<(Fine, bool)>> = Vec::new();
    for cell in &cells {
        let mut hs = Vec::new();
        match *cell {
            Cell::Wall(g, s) => {
                let ml = material_left[&(g, s)];
                match ends(g) {
                    Some((p, e)) if ml => hs.extend([
                        (Fine::H(g, s), true),
                        (Fine::V(e, s), true),
                        (Fine::H(g, s + 1), false),
                        (Fine::V(p, s), false),
                    ]),
                    Some((p, e)) => hs.extend([
                        (Fine::H(g, s), false),
                        (Fine::V(p, s), true),
                        (Fine::H(g, s + 1), true),
                        (Fine::V(e, s), false),
                    ]),
                    None => hs.extend([(Fine::H(g, s), ml), (Fine::H(g, s + 1), !ml)]),
                }
            }
            Cell::Cap(k, up, i) => {
                for cycle in &caps[&(k, up)][i].1 {
                    for &(g, forward) in cycle {
                        hs.push((Fine::H(g, k), forward == up));
                    }
                }
            }
        }
        half.push(hs);
    }
    let mut uses: BTreeMap<Fine, Vec<(usize, bool)>> = BTreeMap::new();
    for (c, hs) in half.iter().enumerate() {
        for &(f, d) in hs {
            uses.entry(f).or_default().push((c, d));
        }
    }
    for u in uses.values() {
        match u.len() {
            2 if u[0].1 != u[1].1 => {}
            4 => return Err(Error::Degenerate("a result touching itself along an edge")),
            _ => return Err(Error::InvalidTopology("an open or misoriented stack")),
        }
    }
    // Faces: walls joined across fine edges where they run on.
    let wall_of = |c: usize| match cells[c] {
        Cell::Wall(g, s) => Some((g, s)),
        Cell::Cap(..) => None,
    };
    // Whether a wall's region-left direction runs counter-clockwise (on a
    // circle) or its direction (on a line).
    let region_left = |g: usize, s: usize| material_left[&(g, s)];
    let mut parent: Vec<usize> = (0..cells.len()).collect();
    let mut interior: BTreeSet<Fine> = BTreeSet::new();
    for (f, u) in &uses {
        let (Some((g1, s1)), Some((g2, s2))) = (wall_of(u[0].0), wall_of(u[1].0)) else {
            continue;
        };
        let joins = match f {
            Fine::H(..) => true,
            Fine::V(..) => {
                let (x, y) = (&views[g1], &views[g2]);
                same_carrier(x, y)
                    && match circle(x) {
                        Some(_) => {
                            (stored_ccw(x) == region_left(g1, s1))
                                == (stored_ccw(y) == region_left(g2, s2))
                        }
                        None => {
                            let d = |v: &PieceView, ml: bool| {
                                let (a, b) = if ml { (v.p, v.e) } else { (v.e, v.p) };
                                [q(b.x) - q(a.x), q(b.y) - q(a.y)]
                            };
                            let (a, b) = (d(x, region_left(g1, s1)), d(y, region_left(g2, s2)));
                            &a[0] * &b[0] + &a[1] * &b[1] > q(0.0)
                        }
                    }
            }
        };
        if joins {
            union(&mut parent, u[0].0, u[1].0);
            interior.insert(*f);
        }
    }
    // Groups in cell order.
    let mut group_of: Vec<usize> = vec![0; cells.len()];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut root_group: BTreeMap<usize, usize> = BTreeMap::new();
    for (c, group) in group_of.iter_mut().enumerate() {
        let r = find(&mut parent, c);
        let gi = *root_group.entry(r).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[gi].push(c);
        *group = gi;
    }
    for (f, u) in &uses {
        if !interior.contains(f) && group_of[u[0].0] == group_of[u[1].0] {
            return Err(Error::Degenerate("a face touching itself along an edge"));
        }
    }

    // Vertices: where the remaining fine edges do not run straight on
    // between the same two faces.
    let fine_ends = |f: Fine| -> Option<(Node, Node)> {
        match f {
            Fine::H(g, k) => ends(g).map(|(p, e)| ((p, k), (e, k))),
            Fine::V(p, s) => Some(((p, s), (p, s + 1))),
        }
    };
    let kept: Vec<Fine> = uses
        .keys()
        .filter(|f| !interior.contains(f))
        .copied()
        .collect();
    let mut incident: BTreeMap<Node, Vec<Fine>> = BTreeMap::new();
    for &f in &kept {
        if let Some((a, b)) = fine_ends(f) {
            incident.entry(a).or_default().push(f);
            incident.entry(b).or_default().push(f);
        }
    }
    let faces_of =
        |f: Fine| -> BTreeSet<usize> { uses[&f].iter().map(|u| group_of[u.0]).collect() };
    let removable = |node: &Node| -> bool {
        let fs = &incident[node];
        if fs.len() != 2 || faces_of(fs[0]) != faces_of(fs[1]) {
            return false;
        }
        match (fs[0], fs[1]) {
            (Fine::V(..), Fine::V(..)) => true,
            (Fine::H(g1, _), Fine::H(g2, _)) => same_carrier(&views[g1], &views[g2]),
            _ => false,
        }
    };
    let mut parts = TopologyParts::default();
    let mut vertex_of: BTreeMap<Node, VertexId> = BTreeMap::new();
    for node in incident.keys() {
        if !removable(node) {
            vertex_of.insert(*node, VertexId(parts.vertices.len()));
            parts.vertices.push(Vertex {
                position: at(node.0, node.1),
                enclosure: None,
            });
        }
    }
    // Edges: chains of fine edges between vertices, and rings.
    let mut edge_of: BTreeMap<Fine, (EdgeId, bool)> = BTreeMap::new();
    let mut chains: Vec<Vec<(Fine, bool)>> = Vec::new();
    let walk = |start: Node, first: Fine, taken: &BTreeMap<Fine, (EdgeId, bool)>| {
        let mut chain = Vec::new();
        let (mut node, mut f) = (start, first);
        loop {
            let (a, b) = fine_ends(f).expect("a fine edge with ends");
            let forward = a == node;
            chain.push((f, forward));
            node = if forward { b } else { a };
            if node == start || !incident.contains_key(&node) || vertex_of.contains_key(&node) {
                break;
            }
            let next = incident[&node]
                .iter()
                .copied()
                .find(|x| *x != f && !taken.contains_key(x));
            match next {
                Some(x) => f = x,
                None => break,
            }
        }
        chain
    };
    for node in vertex_of.keys().copied().collect::<Vec<_>>() {
        for f in incident[&node].clone() {
            if edge_of.contains_key(&f) {
                continue;
            }
            let chain = walk(node, f, &edge_of);
            let id = EdgeId(chains.len());
            for &(x, d) in &chain {
                edge_of.insert(x, (id, d));
            }
            chains.push(chain);
        }
    }
    for &f in &kept {
        if edge_of.contains_key(&f) {
            continue;
        }
        let chain = match fine_ends(f) {
            None => vec![(f, true)],
            Some((a, _)) => walk(a, f, &edge_of),
        };
        let id = EdgeId(chains.len());
        for &(x, d) in &chain {
            edge_of.insert(x, (id, d));
        }
        chains.push(chain);
    }
    // Each chain's curve, oriented upward or counter-clockwise.
    for (i, chain) in chains.iter_mut().enumerate() {
        let (f0, d0) = chain[0];
        let flip = match f0 {
            Fine::V(..) => !d0,
            Fine::H(g, _) => circle(&views[g]).is_some() && stored_ccw(&views[g]) != d0,
        };
        if flip {
            chain.reverse();
            for x in chain.iter_mut() {
                x.1 = !x.1;
            }
            for &(x, d) in chain.iter() {
                edge_of.insert(x, (EdgeId(i), d));
            }
        }
        let node_at = |(f, d): (Fine, bool), start: bool| -> Option<Node> {
            fine_ends(f).map(|(a, b)| if d == start { a } else { b })
        };
        let (first, last) = (chain[0], *chain.last().expect("a chain"));
        let (s, e) = (node_at(first, true), node_at(last, false));
        let (sv, ev) = (
            s.and_then(|x| vertex_of.get(&x).copied()),
            e.and_then(|x| vertex_of.get(&x).copied()),
        );
        let curve = match first.0 {
            Fine::V(..) => Curve3::LineSegment {
                start: parts.vertices[sv.expect("a vertex").0].position,
                end: parts.vertices[ev.expect("a vertex").0].position,
            },
            Fine::H(g, k) => match circle(&views[g]) {
                None => Curve3::LineSegment {
                    start: parts.vertices[sv.expect("a vertex").0].position,
                    end: parts.vertices[ev.expect("a vertex").0].position,
                },
                Some((center, radius)) => {
                    let ring_frame =
                        Frame3::new(frame.point(center, heights[k]), n, frame.x(), tolerance)?;
                    if sv.is_none() {
                        Curve3::Circle {
                            frame: ring_frame,
                            radius,
                        }
                    } else {
                        let mut sweep = 0.0;
                        for &(f, d) in chain.iter() {
                            let Fine::H(h, _) = f else {
                                unreachable!("a horizontal chain")
                            };
                            let (a, b) = if d {
                                (views[h].p, views[h].e)
                            } else {
                                (views[h].e, views[h].p)
                            };
                            sweep += crate::profile::arc_sweep(center, a, b, true);
                        }
                        let p = arranged.position(s.expect("a start").0);
                        Curve3::CircularArc {
                            frame: ring_frame,
                            radius,
                            start_angle: (p.y - center.y).atan2(p.x - center.x),
                            sweep_angle: sweep,
                        }
                    }
                }
            },
        };
        parts.edges.push(Edge {
            start: sv,
            end: ev,
            curve,
            fins: Vec::new(),
        });
    }

    // Faces: surfaces, fins and loops.
    struct Built {
        surface: Surface,
        sense: Orientation,
        loops: Vec<(Vec<Fin>, [i32; 2])>,
    }
    let mut built: Vec<Built> = Vec::new();
    for group in &groups {
        let first = cells[group[0]];
        let low = group
            .iter()
            .map(|&c| match cells[c] {
                Cell::Wall(_, s) => s,
                Cell::Cap(k, ..) => k,
            })
            .min()
            .expect("a cell");
        let (surface, sense, cylinder) = match first {
            Cell::Cap(k, up, _) => {
                let normal = if up { n } else { -n };
                let plane = Frame3::new(
                    frame.point(Point2::default(), heights[k]),
                    normal,
                    frame.x(),
                    tolerance,
                )?;
                (Surface::Plane(plane), Orientation::Forward, None)
            }
            Cell::Wall(g, s) => {
                let v = &views[g];
                let ml = material_left[&(g, s)];
                match circle(v) {
                    None => {
                        let (a, b) = if ml { (v.p, v.e) } else { (v.e, v.p) };
                        let tangent = frame.x() * (b.x - a.x) + frame.y() * (b.y - a.y);
                        let plane = Frame3::new(
                            frame.point(a, heights[low]),
                            tangent.cross(n),
                            tangent,
                            tolerance,
                        )?;
                        (Surface::Plane(plane), Orientation::Forward, None)
                    }
                    Some((center, radius)) => {
                        let axis = Frame3::new(
                            frame.point(center, heights[low]),
                            n,
                            frame.x(),
                            tolerance,
                        )?;
                        // Its normal leaves the material: outward when the
                        // material lies inside, run counter-clockwise.
                        let sense = if stored_ccw(v) == ml {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        (
                            Surface::Cylinder {
                                frame: axis,
                                radius,
                            },
                            sense,
                            Some((center, heights[low])),
                        )
                    }
                }
            }
        };
        // The face's uses of edges.
        let mut used: BTreeMap<EdgeId, Orientation> = BTreeMap::new();
        for &c in group {
            for &(f, d) in &half[c] {
                if interior.contains(&f) {
                    continue;
                }
                let (e, along) = edge_of[&f];
                let o = if d == along {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if used.insert(e, o).is_some_and(|x| x != o) {
                    return Err(Error::Degenerate("a face using an edge both ways"));
                }
            }
        }
        // Loops: each fin's end starts exactly one fin.
        let fin_ends = |e: EdgeId, o: Orientation| {
            let edge = &parts.edges[e.0];
            if o == Orientation::Forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            }
        };
        let mut starts: BTreeMap<VertexId, EdgeId> = BTreeMap::new();
        for (&e, &o) in &used {
            if let (Some(a), _) = fin_ends(e, o) {
                if starts.insert(a, e).is_some() {
                    return Err(Error::Degenerate("a face touching itself at a vertex"));
                }
            }
        }
        let mut done: BTreeSet<EdgeId> = BTreeSet::new();
        let mut loops: Vec<Vec<(EdgeId, Orientation)>> = Vec::new();
        for (&e, &o) in &used {
            if done.contains(&e) {
                continue;
            }
            done.insert(e);
            let mut cycle = vec![(e, o)];
            let (start, mut to) = fin_ends(e, o);
            while to != start {
                let next = to
                    .and_then(|v| starts.get(&v))
                    .copied()
                    .ok_or(Error::InvalidTopology("an open face loop"))?;
                if !done.insert(next) {
                    return Err(Error::InvalidTopology("a face loop does not close"));
                }
                cycle.push((next, used[&next]));
                to = fin_ends(next, used[&next]).1;
            }
            loops.push(cycle);
        }
        // Pcurves.
        let mut out: Vec<(Vec<Fin>, [i32; 2], f64)> = Vec::new();
        for cycle in loops {
            let mut fins = Vec::new();
            match (&surface, cylinder) {
                (Surface::Plane(plane), _) => {
                    let mut twice = 0.0;
                    for &(e, o) in &cycle {
                        let pcurve = plane_pcurve(&parts.edges[e.0].curve, o, *plane);
                        for i in 0..16 {
                            let (a, b) = (
                                pcurve.point(f64::from(i) / 16.0),
                                pcurve.point(f64::from(i + 1) / 16.0),
                            );
                            twice += a.x * b.y - b.x * a.y;
                        }
                        fins.push(Fin {
                            edge: e,
                            sense: o,
                            pcurve,
                            enclosure: None,
                        });
                    }
                    out.push((fins, [0, 0], twice));
                }
                (_, Some((center, base))) => {
                    // In the cover: u continues along the loop.
                    let angle = |p: Point2| (p.y - center.y).atan2(p.x - center.x);
                    let first = cycle[0];
                    let mut u = match fin_ends(first.0, first.1).0 {
                        Some(v) => {
                            let [x, y, _] = frame.coordinates(parts.vertices[v.0].position);
                            angle(Point2::new(x, y))
                        }
                        None if first.1 == Orientation::Forward => 0.0,
                        None => TAU,
                    };
                    let u0 = u;
                    for &(e, o) in &cycle {
                        let edge = &parts.edges[e.0];
                        let sign = if o == Orientation::Forward { 1.0 } else { -1.0 };
                        let (from, to) = match &edge.curve {
                            Curve3::LineSegment { start, end } => {
                                let (a, b) = if o == Orientation::Forward {
                                    (*start, *end)
                                } else {
                                    (*end, *start)
                                };
                                let (va, vb) = (
                                    frame.coordinates(a)[2] - base,
                                    frame.coordinates(b)[2] - base,
                                );
                                (Point2::new(u, va), Point2::new(u, vb))
                            }
                            Curve3::CircularArc {
                                frame: arc,
                                sweep_angle,
                                ..
                            } => {
                                let v = frame.coordinates(arc.origin())[2] - base;
                                let next = u + sign * sweep_angle;
                                let pair = (Point2::new(u, v), Point2::new(next, v));
                                u = next;
                                pair
                            }
                            Curve3::Circle { frame: ring, .. } => {
                                let v = frame.coordinates(ring.origin())[2] - base;
                                let next = u + sign * TAU;
                                let pair = (Point2::new(u, v), Point2::new(next, v));
                                u = next;
                                pair
                            }
                            _ => return Err(Error::InvalidTopology("a stack wall's edge")),
                        };
                        fins.push(Fin {
                            edge: e,
                            sense: o,
                            pcurve: Curve2::LineSegment {
                                start: from,
                                end: to,
                            },
                            enclosure: None,
                        });
                    }
                    let turns = ((u - u0) / TAU).round() as i32;
                    // Its area in the cover, against the face's sense (a
                    // loop winding round the axis first).
                    let mut twice = 0.0;
                    for fin in &fins {
                        if let Curve2::LineSegment { start: a, end: b } = fin.pcurve {
                            twice += a.x * b.y - b.x * a.y;
                        }
                    }
                    let key = if turns != 0 {
                        f64::INFINITY
                    } else if sense == Orientation::Forward {
                        twice
                    } else {
                        -twice
                    };
                    out.push((fins, [turns, 0], key));
                }
                _ => unreachable!("a plane or a cylinder"),
            }
        }
        // The outer loop (counter-clockwise about the face's normal) comes
        // first.
        out.sort_by(|a, b| b.2.total_cmp(&a.2));
        // On a cylinder every loop lies on the first's sheet of the cover:
        // one not winding round starts within a turn of its lowest angle.
        if cylinder.is_some() && out.len() > 1 {
            let lowest = out[0]
                .0
                .iter()
                .filter_map(|f| match f.pcurve {
                    Curve2::LineSegment { start, end } => Some(start.x.min(end.x)),
                    _ => None,
                })
                .fold(f64::INFINITY, f64::min);
            for (fins, winding, _) in out.iter_mut().skip(1) {
                let Some(Curve2::LineSegment { start, .. }) =
                    fins.first().map(|f| f.pcurve.clone())
                else {
                    continue;
                };
                let turns = ((start.x - lowest) / TAU).floor();
                if winding[0] != 0 || turns == 0.0 {
                    continue;
                }
                for fin in fins.iter_mut() {
                    if let Curve2::LineSegment { start, end } = &mut fin.pcurve {
                        start.x -= turns * TAU;
                        end.x -= turns * TAU;
                    }
                }
            }
        }
        built.push(Built {
            surface,
            sense,
            loops: out.into_iter().map(|(f, w, _)| (f, w)).collect(),
        });
    }

    // Components: faces joined by edges.
    let mut face_parent: Vec<usize> = (0..groups.len()).collect();
    for &f in &kept {
        let u = &uses[&f];
        union(&mut face_parent, group_of[u[0].0], group_of[u[1].0]);
    }
    let mut members: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for gi in 0..groups.len() {
        let r = find(&mut face_parent, gi);
        members.entry(r).or_default().push(gi);
    }
    // A shell whose lowest face faces up bounds a cavity.
    let lowest_up = |faces: &[usize]| -> bool {
        faces
            .iter()
            .filter_map(|&gi| match cells[groups[gi][0]] {
                Cell::Cap(k, up, _) => Some((k, up)),
                Cell::Wall(..) => None,
            })
            .min_by_key(|(k, up)| (*k, *up))
            .is_some_and(|(_, up)| up)
    };
    let (mut outers, mut cavities): (Vec<Vec<usize>>, Vec<Vec<usize>>) = (Vec::new(), Vec::new());
    for faces in members.into_values() {
        if lowest_up(&faces) {
            cavities.push(faces);
        } else {
            outers.push(faces);
        }
    }
    if !cavities.is_empty() && outers.len() != 1 {
        return Err(Error::OutOfDomain("a stack's cavity among several solids"));
    }

    // Plans' inputs: every stored vertex's point, and each segment's ends.
    let profiles = [(Operand::A, &stack.a), (Operand::B, &stack.b)];
    let mut owners: BTreeMap<PId, Vec<(Operand, usize, usize)>> = BTreeMap::new();
    let mut seg_ends: BTreeMap<(Operand, usize, usize), [PId; 2]> = BTreeMap::new();
    for (o, profile) in profiles {
        for (b, boundary) in profile.boundaries().enumerate() {
            let count = boundary
                .polygon_vertices()
                .or(boundary.path_geometry().map(|g| g.0))
                .map_or(0, |v| v.len());
            for j in 0..count {
                let p = arranged.canon(PId::Vertex(o, b, j));
                owners.entry(p).or_default().push((o, b, j));
                seg_ends.insert(
                    (o, b, j),
                    [p, arranged.canon(PId::Vertex(o, b, (j + 1) % count))],
                );
            }
        }
    }
    // Whether a stored segment ends at a point (a whole circle has no
    // ends).
    let ends_at = |o: Operand, b: usize, j: usize, p: PId| {
        seg_ends.get(&(o, b, j)).is_some_and(|e| e.contains(&p))
    };
    let range = |o: Operand| if o == Operand::A { stack.ha } else { stack.hb };
    let end_at = |o: Operand, h: f64| -> Option<At> {
        let r = range(o);
        if h == r[0] {
            Some(At::Low)
        } else if h == r[1] {
            Some(At::High)
        } else {
            None
        }
    };
    let within = |o: Operand, h: f64| {
        let r = range(o);
        r[0] < h && h < r[1]
    };
    let spans = |o: Operand, s: usize| {
        let r = range(o);
        r[0] <= heights[s] && heights[s + 1] <= r[1]
    };
    let tool = |o: Operand| stack.op == Op2::Cut && o == Operand::B;
    let keep = |o: Operand, key: Key, continues: &mut Vec<Key>, touches: &mut Vec<Key>| {
        if tool(o) {
            touches.push(key);
        } else {
            continues.push(key);
        }
    };
    // The stored segments a piece lies on: its own and its partner's.
    let segs_of = |g: usize| -> Vec<(Operand, usize, usize)> {
        let v = &views[g];
        let mut out = vec![(v.op, v.b, v.j)];
        if let crate::profile::boolean::Class::Shared { partner, .. } = v.class {
            let w = &views[partner];
            out.push((w.op, w.b, w.j));
        }
        out
    };
    // Every stored segment through or ending at each point.
    let mut through: BTreeMap<PId, BTreeSet<(Operand, usize, usize)>> = BTreeMap::new();
    for v in &views {
        for p in [v.from, v.to].into_iter().flatten() {
            through.entry(p).or_default().insert((v.op, v.b, v.j));
        }
    }
    // Cap faces' plans first: edges and vertices inside a wall touch them.
    let other_profile = |o: Operand| if o == Operand::A { &stack.a } else { &stack.b };
    let mut cap_plan: BTreeMap<usize, (Vec<Key>, Vec<Key>)> = BTreeMap::new();
    for (gi, group) in groups.iter().enumerate() {
        let Cell::Cap(k, up, i) = cells[group[0]] else {
            continue;
        };
        let (mut continues, mut touches) = (Vec::new(), Vec::new());
        for o in [Operand::A, Operand::B] {
            let Some(e) = end_at(o, heights[k]) else {
                continue;
            };
            let facing_up = e == At::High;
            let overlaps = !boolean(&caps[&(k, up)][i].0, other_profile(o), Op2::Common)?
                .pieces
                .is_empty();
            if !overlaps {
                continue;
            }
            let key = (o, e, What::Cap, 0, 0);
            if facing_up == up && !tool(o) {
                continues.push(key);
            } else if tool(o) && facing_up != up {
                touches.push(key);
            }
        }
        cap_plan.insert(gi, (continues, touches));
    }
    // The caps an edge's or vertex's faces lie on.
    let caps_touched = |faces: &BTreeSet<usize>| -> Vec<Key> {
        faces
            .iter()
            .filter_map(|gi| cap_plan.get(gi))
            .flat_map(|(c, t)| c.iter().chain(t).copied())
            .collect()
    };

    let mut out = Vec::new();
    for (index, faces) in outers.iter().enumerate() {
        let cavity: Vec<usize> = if index == 0 {
            cavities.iter().flatten().copied().collect()
        } else {
            Vec::new()
        };
        let mut all: Vec<usize> = faces.clone();
        all.extend(&cavity);
        all.sort_unstable();
        let face_id: BTreeMap<usize, FaceId> = all
            .iter()
            .enumerate()
            .map(|(i, gi)| (*gi, FaceId(i)))
            .collect();
        // Its edges and vertices, renumbered in order.
        let mut edges: BTreeSet<EdgeId> = BTreeSet::new();
        for &gi in &all {
            for (fins, _) in &built[gi].loops {
                edges.extend(fins.iter().map(|f| f.edge));
            }
        }
        let edge_id: BTreeMap<EdgeId, EdgeId> = edges
            .iter()
            .enumerate()
            .map(|(i, e)| (*e, EdgeId(i)))
            .collect();
        let mut vertices: BTreeSet<VertexId> = BTreeSet::new();
        for e in &edges {
            let edge = &parts.edges[e.0];
            vertices.extend(edge.start);
            vertices.extend(edge.end);
        }
        let vertex_id: BTreeMap<VertexId, VertexId> = vertices
            .iter()
            .enumerate()
            .map(|(i, v)| (*v, VertexId(i)))
            .collect();
        let mut p = TopologyParts::default();
        for v in &vertices {
            p.vertices.push(parts.vertices[v.0].clone());
        }
        for e in &edges {
            let edge = &parts.edges[e.0];
            p.edges.push(Edge {
                start: edge.start.map(|v| vertex_id[&v]),
                end: edge.end.map(|v| vertex_id[&v]),
                curve: edge.curve.clone(),
                fins: Vec::new(),
            });
        }
        let has_cavity = !cavity.is_empty();
        for &gi in &all {
            let inner = cavity.contains(&gi);
            let mut loop_ids = Vec::new();
            for (fins, winding) in &built[gi].loops {
                let mut ids = Vec::new();
                for fin in fins {
                    let mut fin = fin.clone();
                    fin.edge = edge_id[&fin.edge];
                    let id = FinId(p.fins.len());
                    p.edges[fin.edge.0].fins.push(id);
                    p.fins.push(fin);
                    ids.push(id);
                }
                loop_ids.push(LoopId(p.loops.len()));
                p.loops.push(Loop::Edges {
                    fins: ids,
                    winding: *winding,
                });
            }
            p.faces.push(Face {
                surface: built[gi].surface.clone(),
                sense: built[gi].sense,
                loops: loop_ids,
                front: if inner { ShellId(2) } else { ShellId(0) },
                back: if inner { ShellId(3) } else { ShellId(1) },
                enclosure: None,
            });
        }
        let sides = |inner: bool, side: Side| -> Vec<(FaceId, Side)> {
            all.iter()
                .filter(|gi| cavity.contains(gi) == inner)
                .map(|gi| (face_id[gi], side))
                .collect()
        };
        p.shells = vec![
            Shell {
                region: RegionId(1),
                sides: sides(false, Side::Front),
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: sides(false, Side::Back),
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
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
            p.shells.push(Shell {
                region: RegionId(1),
                sides: sides(true, Side::Front),
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            });
            p.shells.push(Shell {
                region: RegionId(2),
                sides: sides(true, Side::Back),
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            });
            p.regions[1].shells.push(ShellId(2));
            p.regions.push(Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(3)],
            });
        }

        // Plans.
        let mut plans: Vec<KeyPlan> = Vec::new();
        let finish = |mut c: Vec<Key>, mut t: Vec<Key>| {
            c.sort();
            c.dedup();
            t.sort();
            t.dedup();
            t.retain(|x| !c.contains(x));
            (c, t)
        };
        // Faces.
        let mut members: BTreeSet<Operand> = BTreeSet::new();
        for &gi in &all {
            let (mut continues, mut touches) = (Vec::new(), Vec::new());
            let role = match cells[groups[gi][0]] {
                Cell::Cap(_, up, _) => {
                    let (c, t) = &cap_plan[&gi];
                    continues.extend(c.iter().copied());
                    touches.extend(t.iter().copied());
                    if up {
                        Role::EndCap
                    } else {
                        Role::StartCap
                    }
                }
                Cell::Wall(..) => {
                    for &c in &groups[gi] {
                        let (g, s) = wall_of(c).expect("a wall");
                        let ml = material_left[&(g, s)];
                        let v = &views[g];
                        for (o, b, j) in segs_of(g) {
                            if !spans(o, s) {
                                continue;
                            }
                            let own = if o == Operand::A {
                                v.left[0]
                            } else {
                                v.left[1]
                            };
                            let key = (o, At::Swept, What::Wall, b, j);
                            if own == ml && !tool(o) {
                                continues.push(key);
                            } else {
                                touches.push(key);
                            }
                        }
                    }
                    Role::Wall
                }
            };
            members.extend(continues.iter().map(|k| k.0));
            let (c, t) = finish(continues, touches);
            plans.push((Slot::Face(face_id[&gi]), c, t, EntityKind::Face, role));
        }
        // Edges.
        let mut edge_role: BTreeMap<EdgeId, Role> = BTreeMap::new();
        for e in &edges {
            let chain = &chains[e.0];
            let faces: BTreeSet<usize> = chain.iter().flat_map(|(f, _)| faces_of(*f)).collect();
            let (mut continues, mut touches) = (Vec::new(), Vec::new());
            let role = match chain[0].0 {
                Fine::V(..) => {
                    for &(f, _) in chain {
                        let Fine::V(v, s) = f else { continue };
                        for &(o, b, j) in owners.get(&v).into_iter().flatten() {
                            if spans(o, s) {
                                keep(
                                    o,
                                    (o, At::Swept, What::Vertical, b, j),
                                    &mut continues,
                                    &mut touches,
                                );
                            }
                        }
                        for &(o, b, j) in through.get(&v).into_iter().flatten() {
                            if spans(o, s) && !ends_at(o, b, j, v) {
                                touches.push((o, At::Swept, What::Wall, b, j));
                            }
                        }
                    }
                    Role::Vertical
                }
                Fine::H(_, k) => {
                    for &(f, _) in chain {
                        let Fine::H(g, _) = f else { continue };
                        for (o, b, j) in segs_of(g) {
                            match end_at(o, heights[k]) {
                                Some(e) => {
                                    keep(o, (o, e, What::Edge, b, j), &mut continues, &mut touches)
                                }
                                None if within(o, heights[k]) => {
                                    touches.push((o, At::Swept, What::Wall, b, j));
                                    touches.extend(caps_touched(&faces));
                                }
                                None => {}
                            }
                        }
                    }
                    let up = faces
                        .iter()
                        .any(|gi| matches!(cells[groups[*gi][0]], Cell::Cap(_, true, _)));
                    if up {
                        Role::TopEdge
                    } else {
                        Role::BottomEdge
                    }
                }
            };
            edge_role.insert(*e, role);
            let (c, t) = finish(continues, touches);
            plans.push((Slot::Edge(edge_id[e]), c, t, EntityKind::Edge, role));
        }
        // Vertices.
        for v in &vertices {
            let node = *vertex_of
                .iter()
                .find(|(_, x)| *x == v)
                .map(|(n, _)| n)
                .expect("a vertex's node");
            let (point, k) = node;
            let h = heights[k];
            let faces: BTreeSet<usize> =
                incident[&node].iter().flat_map(|f| faces_of(*f)).collect();
            let (mut continues, mut touches) = (Vec::new(), Vec::new());
            for &(o, b, j) in owners.get(&point).into_iter().flatten() {
                match end_at(o, h) {
                    Some(e) => keep(o, (o, e, What::Vertex, b, j), &mut continues, &mut touches),
                    None if within(o, h) => {
                        touches.push((o, At::Swept, What::Vertical, b, j));
                        touches.extend(caps_touched(&faces));
                    }
                    None => {}
                }
            }
            for &(o, b, j) in through.get(&point).into_iter().flatten() {
                if ends_at(o, b, j, point) {
                    continue;
                }
                match end_at(o, h) {
                    Some(e) => touches.push((o, e, What::Edge, b, j)),
                    None if within(o, h) => touches.push((o, At::Swept, What::Wall, b, j)),
                    None => {}
                }
            }
            let role = if incident[&node].iter().any(|f| {
                matches!(f, Fine::H(..))
                    && edge_of
                        .get(f)
                        .and_then(|(e, _)| edge_role.get(e))
                        .is_some_and(|r| *r == Role::TopEdge)
            }) {
                Role::TopVertex
            } else {
                Role::BottomVertex
            };
            let (c, t) = finish(continues, touches);
            plans.push((Slot::Vertex(vertex_id[v]), c, t, EntityKind::Vertex, role));
        }
        // Regions: the solid continues the inputs whose material it holds;
        // a cavity is new, where the tool was.
        if stack.op == Op2::Cut {
            members = BTreeSet::from([Operand::A]);
        }
        let regions: Vec<Key> = members
            .iter()
            .map(|o| (*o, At::Swept, What::Region, 0, 0))
            .collect();
        plans.push((
            Slot::Region(RegionId(1)),
            regions,
            Vec::new(),
            EntityKind::Region,
            Role::Region,
        ));
        if has_cavity {
            plans.push((
                Slot::Region(RegionId(2)),
                Vec::new(),
                vec![(Operand::B, At::Swept, What::Region, 0, 0)],
                EntityKind::Region,
                Role::Region,
            ));
        }
        // Height range.
        let levels: Vec<usize> = all
            .iter()
            .filter_map(|gi| match cells[groups[*gi][0]] {
                Cell::Cap(k, ..) => Some(k),
                Cell::Wall(..) => None,
            })
            .collect();
        let (lo, hi) = (
            *levels.iter().min().expect("a cap"),
            *levels.iter().max().expect("a cap"),
        );
        out.push(Component {
            parts: p,
            plans,
            heights: [heights[lo], heights[hi]],
        });
    }
    Ok(out)
}
