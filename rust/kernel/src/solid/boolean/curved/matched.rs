//! S9e.2: a given result whose stored topology comes from another assembly
//! (a stack's slabs, S9b.1's fragments) or holds one solid of several
//! (REVIEW_NOTES.md, S9e.2's refined decisions).
//!
//! The re-run's assembly is matched to the stored topology geometrically:
//! each stored vertex the one re-run vertex within the resolution of it,
//! one to one and onto; each stored edge the re-run edge between the
//! matched ends whose points at a quarter, a half and three quarters of its
//! range lie within the resolution of the stored edge's at those fractions
//! (either direction; a ring's points within the resolution of the stored
//! ring's curve); each stored face the re-run face bounded by the matched
//! edges on a surface of the same kind. A result of several solids is given
//! as its whole construction, and the second arrangement's pieces are
//! sorted by solid once it is built, by adjacency (`keep_solid`): pieces of
//! the other solids' faces are dropped, and the other input's pieces inside
//! another solid are outside the given one.
use super::graph::{Arr, CurveRef};
use crate::topology::{Curve3, Loop, Topology, TopologyParts, VertexId};
use crate::{Error, Point3, Result};
use std::collections::{BTreeMap, BTreeSet};

/// A re-run solid's slots matched to the stored ones: each re-run face,
/// edge and vertex slot's stored slot.
#[derive(Debug, Clone)]
pub(crate) struct Match {
    pub(crate) faces: Vec<usize>,
    pub(crate) edges: Vec<usize>,
    pub(crate) vertices: Vec<usize>,
}

fn near(a: Point3, b: Point3, tol: f64) -> bool {
    let (a, b) = (a.to_array(), b.to_array());
    (0..3).all(|k| (a[k] - b[k]).abs() <= tol)
}

/// A point's distance from a ring's closed curve (a circle or a whole
/// ellipse): its offset from the curve's plane and from the curve in it,
/// the ellipse's point at the point's eccentric angle (an upper bound).
fn ring_distance(curve: &Curve3, p: Point3) -> Option<f64> {
    match curve {
        Curve3::Circle { frame, radius } | Curve3::CircularArc { frame, radius, .. } => {
            let [x, y, z] = frame.coordinates(p);
            Some(((x.hypot(y) - radius).powi(2) + z * z).sqrt())
        }
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            ..
        } => {
            let [x, y, z] = frame.coordinates(p);
            let t = (y / minor).atan2(x / major);
            let (s, c) = t.sin_cos();
            Some(((x - major * c).powi(2) + (y - minor * s).powi(2) + z * z).sqrt())
        }
        _ => None,
    }
}

/// A face's bounding edges and vertex loops.
fn bounds_of(
    loops: &[crate::topology::LoopId],
    all: &[Loop],
    fins: &[crate::topology::Fin],
) -> (Vec<usize>, Vec<usize>) {
    let (mut edges, mut verts) = (Vec::new(), Vec::new());
    for l in loops {
        match &all[l.index()] {
            Loop::Edges { fins: fs, .. } => {
                edges.extend(fs.iter().map(|f| fins[f.index()].edge.index()));
            }
            Loop::Vertex(v) => verts.push(v.index()),
        }
    }
    edges.sort_unstable();
    edges.dedup();
    verts.sort_unstable();
    (edges, verts)
}

/// The re-run solid `parts` matched to the stored topology `t` within the
/// resolution `tol`, or none where they differ.
pub(crate) fn matched(parts: &TopologyParts, t: &Topology, tol: f64) -> Option<Match> {
    if parts.vertices.len() != t.vertices().len()
        || parts.edges.len() != t.edges().len()
        || parts.faces.len() != t.faces().len()
        || parts.shells.len() != t.shells().len()
    {
        return None;
    }
    // Vertices: one stored vertex within the resolution, one to one.
    let mut vertices = Vec::with_capacity(parts.vertices.len());
    let mut used = BTreeSet::new();
    for v in &parts.vertices {
        let mut hits = t
            .vertices()
            .iter()
            .enumerate()
            .filter(|(_, s)| near(v.position, s.position, tol))
            .map(|(j, _)| j);
        let (Some(j), None) = (hits.next(), hits.next()) else {
            return None;
        };
        if !used.insert(j) {
            return None;
        }
        vertices.push(j);
    }
    // Edges: the stored edge between the matched ends through the re-run
    // edge's points.
    let map = |v: Option<VertexId>| v.map(|v| vertices[v.index()]);
    let mut edges = Vec::with_capacity(parts.edges.len());
    let mut used = BTreeSet::new();
    let fractions = [0.25, 0.5, 0.75];
    for e in &parts.edges {
        let (s, f) = (map(e.start), map(e.end));
        let fits = |x: &crate::topology::Edge| -> bool {
            let (xs, xf) = (x.start.map(VertexId::index), x.end.map(VertexId::index));
            if x.is_ring() || e.is_ring() {
                return x.is_ring()
                    && e.is_ring()
                    && [0.0, 0.25, 0.5, 0.75].iter().all(|&k| {
                        ring_distance(&x.curve, e.curve.point(k)).is_some_and(|d| d <= tol)
                    });
            }
            let along = |rev: bool| {
                fractions.iter().all(|&k| {
                    let q = x.curve.point(if rev { 1.0 - k } else { k });
                    near(e.curve.point(k), q, tol)
                })
            };
            (xs == s && xf == f && along(false)) || (xs == f && xf == s && along(true))
        };
        let mut hits = t
            .edges()
            .iter()
            .enumerate()
            .filter(|(_, x)| fits(x))
            .map(|(j, _)| j);
        let (Some(j), None) = (hits.next(), hits.next()) else {
            return None;
        };
        if !used.insert(j) {
            return None;
        }
        edges.push(j);
    }
    // Faces: the stored face bounded by the matched edges on a surface of
    // the same kind.
    let mut faces = Vec::with_capacity(parts.faces.len());
    let mut used = BTreeSet::new();
    let stored: Vec<(Vec<usize>, Vec<usize>)> = t
        .faces()
        .iter()
        .map(|f| bounds_of(&f.loops, t.loops(), t.fins()))
        .collect();
    for f in &parts.faces {
        let (es, vs) = bounds_of(&f.loops, &parts.loops, &parts.fins);
        let mut es: Vec<usize> = es.iter().map(|&e| edges[e]).collect();
        es.sort_unstable();
        let mut vs: Vec<usize> = vs.iter().map(|&v| vertices[v]).collect();
        vs.sort_unstable();
        let kind = std::mem::discriminant(&f.surface);
        let mut hits = t.faces().iter().enumerate().filter(|(j, x)| {
            std::mem::discriminant(&x.surface) == kind && stored[*j] == (es.clone(), vs.clone())
        });
        let (Some((j, _)), None) = (hits.next(), hits.next()) else {
            return None;
        };
        if !used.insert(j) {
            return None;
        }
        faces.push(j);
    }
    Some(Match {
        faces,
        edges,
        vertices,
    })
}

/// The solid of each piece of a given model's faces: the group of pieces of
/// its model face joined across the edges they share, named by a given edge
/// on the group (every group is bounded by given edges of one solid).
fn given_pieces(arr: &Arr, o: usize, solids: &[usize]) -> Result<BTreeMap<usize, usize>> {
    let mut out = BTreeMap::new();
    let mut by_face: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, p) in arr.pieces.iter().enumerate() {
        if p.op == o {
            by_face.entry(p.face).or_default().push(i);
        }
    }
    for ps in by_face.values() {
        let mut parent: Vec<usize> = (0..ps.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (k, &pi) in ps.iter().enumerate() {
            for &(g, _) in arr.pieces[pi].loops.iter().flatten() {
                users.entry(g).or_default().push(k);
            }
        }
        for us in users.values() {
            for w in us.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                parent[a.max(b)] = a.min(b);
            }
        }
        let mut named: BTreeMap<usize, usize> = BTreeMap::new();
        for (k, &pi) in ps.iter().enumerate() {
            for &(g, _) in arr.pieces[pi].loops.iter().flatten() {
                if let CurveRef::Edge(eo, ei) = arr.edges[g].curve {
                    if eo == o {
                        let r = find(&mut parent, k);
                        if *named.entry(r).or_insert(solids[ei]) != solids[ei] {
                            return Err(Error::Degenerate("solids of a given result touching"));
                        }
                    }
                }
            }
        }
        for (k, &pi) in ps.iter().enumerate() {
            let r = find(&mut parent, k);
            let solid = named.get(&r).ok_or(Error::ComputationLimit(
                "a given result's face part without its edges",
            ))?;
            out.insert(pi, *solid);
        }
    }
    Ok(out)
}

/// S9e.2: the second arrangement of a given result of several solids (its
/// whole construction) reduced to the given solid: each piece of the
/// given model's faces sorted by solid (`given_pieces`), each piece of the
/// other input with a side inside the construction by the given pieces
/// sharing an edge with it, else by its neighbours across the other
/// input's own edges inside the construction; the other solids' pieces are
/// dropped and the other input's pieces inside them are outside (both
/// sides).
pub(super) fn keep_solid(arr: &mut Arr) -> Result<()> {
    for o in 0..2 {
        let Some(g) = arr.models[o].given.as_ref() else {
            continue;
        };
        let Some((solids, given)) = g.several() else {
            continue;
        };
        let solids = solids.to_vec();
        let mine = given_pieces(arr, o, &solids)?;
        // The given pieces using each arrangement edge.
        let mut on_edge: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for (&pi, &s) in &mine {
            for &(e, _) in arr.pieces[pi].loops.iter().flatten() {
                on_edge.entry(e).or_default().insert(s);
            }
        }
        let other = 1 - o;
        let inside: Vec<usize> = (0..arr.pieces.len())
            .filter(|&i| arr.pieces[i].op == other && arr.pieces[i].sides != (false, false))
            .collect();
        let mut parent: Vec<usize> = (0..inside.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        // Neighbours across the other input's own edges, both inside.
        let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (k, &pi) in inside.iter().enumerate() {
            for &(e, _) in arr.pieces[pi].loops.iter().flatten() {
                if matches!(arr.edges[e].curve, CurveRef::Edge(eo, _) if eo == other) {
                    users.entry(e).or_default().push(k);
                }
            }
        }
        for us in users.values() {
            for w in us.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                parent[a.max(b)] = a.min(b);
            }
        }
        let mut named: BTreeMap<usize, usize> = BTreeMap::new();
        for (k, &pi) in inside.iter().enumerate() {
            for &(e, _) in arr.pieces[pi].loops.iter().flatten() {
                for &s in on_edge.get(&e).into_iter().flatten() {
                    let r = find(&mut parent, k);
                    if *named.entry(r).or_insert(s) != s {
                        return Err(Error::Degenerate("solids of a given result touching"));
                    }
                }
            }
        }
        let mut outside: BTreeSet<usize> = BTreeSet::new();
        for (k, &pi) in inside.iter().enumerate() {
            let r = find(&mut parent, k);
            let s = named.get(&r).ok_or(Error::ComputationLimit(
                "a solid inside one of a given result's several solids",
            ))?;
            if *s != given {
                outside.insert(pi);
            }
        }
        for pi in outside {
            arr.pieces[pi].sides = (false, false);
        }
        let keep: Vec<bool> = (0..arr.pieces.len())
            .map(|i| arr.pieces[i].op != o || mine.get(&i) == Some(&given))
            .collect();
        let mut k = 0;
        arr.pieces.retain(|_| {
            k += 1;
            keep[k - 1]
        });
    }
    Ok(())
}
