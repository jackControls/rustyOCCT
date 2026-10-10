//! S9e.4b.4c.1: imported polyhedra in the curved engine (REVIEW_NOTES.md,
//! "S9e.4b.4c refined").
//!
//! An imported polyhedron (S9e.4b.2) has no construction: it is decided on
//! its stored vertices, each stored face the triangles of its polygon of
//! stored vertices (S9b.2's stored model, each triangle exactly planar). As
//! a leaf of the curved engine each triangle is a model face on its exact
//! plane (`Facet`), its edges the stored edges and the triangles' diagonals
//! (lines between stored vertices; a diagonal no stored edge, its id none),
//! its vertices the stored ones; a point's membership is the parity of an
//! exact ray's crossings with the triangles (`Mesh::member`), a point on the
//! surface pushed along directions decided by the triangles at it: each one's
//! wedge there (its plane through the point bounded by the rays of its edges
//! at it) crossed or not by the ray from the pushed point, the others by the
//! ray from the point itself. Every triangle of a stored face names its
//! pieces by that face, so the result joins them back into one face across
//! the diagonals (`assemble`'s pieces of one input face).
use super::model::*;
use super::num::*;
use crate::identity::Role;
use crate::profile::boolean::Operand;
use crate::solid::Solid;
use crate::topology::{EdgeId, FaceId, Orientation, RegionId, Slot, VertexId};
use crate::{Error, Frame3, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A triangle of the mesh: its corners counter-clockwise about its outward
/// normal `n` (twice its vector area), each edge's inward normal in its
/// plane (`n x (c[k+1] - c[k])`), and its binary64 box widened past rounding.
#[derive(Debug, Clone)]
struct Tri {
    c: [V; 3],
    n: V,
    w: [V; 3],
    lo: [f64; 3],
    hi: [f64; 3],
}

/// Where a point lies on a closed triangle: inside it, inside its edge `k`
/// (from corner `k` to `k + 1`) or at its corner `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum At {
    Inside,
    Edge(usize),
    Corner(usize),
}

/// An imported polyhedron's closed surface of triangles (its shells, a
/// cavity's among them), each a model face.
#[derive(Debug, Clone)]
pub(super) struct Mesh {
    tris: Vec<Tri>,
    /// Its shells (the stored topology's: an outer one and its cavities').
    shells: usize,
}

/// Ray directions with every coordinate positive (the boxes' slabs).
const RAYS: [[i64; 3]; 6] = [
    [7, 3, 5],
    [2, 11, 13],
    [17, 5, 3],
    [3, 19, 7],
    [23, 29, 31],
    [37, 13, 41],
];

fn ri(n: i64) -> R {
    int(n)
}

impl Tri {
    fn new(c: [V; 3]) -> Self {
        let n = cross(&sub(&c[1], &c[0]), &sub(&c[2], &c[0]));
        let w = [0, 1, 2].map(|k| cross(&n, &sub(&c[(k + 1) % 3], &c[k])));
        let f = |x: &R| crate::solid::split::rational_f64(x);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &c {
            for i in 0..3 {
                lo[i] = lo[i].min(f(&p[i]));
                hi[i] = hi[i].max(f(&p[i]));
            }
        }
        for i in 0..3 {
            let m = 1e-9 * (1.0 + lo[i].abs().max(hi[i].abs()));
            lo[i] -= m;
            hi[i] += m;
        }
        Self { c, n, w, lo, hi }
    }

    /// The plane's value at a point (positive in front).
    fn plane(&self, p: &QV) -> Qd {
        qdot(&qsub(p, &qv(&self.c[0])), &self.n)
    }

    /// The point's side of edge `k`'s line in the plane (positive inside).
    fn side(&self, k: usize, p: &QV) -> Qd {
        qdot(&qsub(p, &qv(&self.c[k])), &self.w[k])
    }

    /// Whether a binary64 point lies within the widened box.
    fn near(&self, x: &[f64; 3]) -> bool {
        (0..3).all(|i| {
            let m = 1e-9 * (1.0 + x[i].abs());
            x[i] >= self.lo[i] - m && x[i] <= self.hi[i] + m
        })
    }

    /// Where a point of its plane lies on the closed triangle, if on it.
    fn at(&self, p: &QV) -> Option<At> {
        let s = [0, 1, 2].map(|k| self.side(k, p).sign());
        if s.contains(&Ordering::Less) {
            return None;
        }
        let zeros: Vec<usize> = (0..3).filter(|&k| s[k] == Ordering::Equal).collect();
        Some(match zeros.as_slice() {
            [] => At::Inside,
            [k] => At::Edge(*k),
            // On edges `k` and `k + 1`: their common corner `k + 1`.
            [a, b] => {
                if (a + 1) % 3 == *b {
                    At::Corner(*b)
                } else {
                    At::Corner(*a)
                }
            }
            _ => return None,
        })
    }

    /// The wedge's bounding normals in its plane at a point `at` it: none
    /// inside, the edge's inward normal on an edge, both edges' at a
    /// corner.
    fn bounds(&self, at: At) -> Vec<&V> {
        match at {
            At::Inside => Vec::new(),
            At::Edge(k) => vec![&self.w[k]],
            At::Corner(k) => vec![&self.w[k], &self.w[(k + 2) % 3]],
        }
    }
}

/// The sign of a linear form along symbolic pushes `d1 + e d2 + ...`: its
/// first nonzero value's.
fn lex(dirs: &[QV], g: &V) -> Ordering {
    for d in dirs {
        let s = qdot(d, g).sign();
        if s != Ordering::Equal {
            return s;
        }
    }
    Ordering::Equal
}

fn times(a: Ordering, b: Ordering) -> Ordering {
    match (a, b) {
        (Ordering::Equal, _) | (_, Ordering::Equal) => Ordering::Equal,
        (x, y) if x == y => Ordering::Greater,
        _ => Ordering::Less,
    }
}

impl Mesh {
    /// Its shells.
    pub(super) fn shells(&self) -> usize {
        self.shells
    }

    /// The triangles holding a point on their closed region, with where.
    fn holding(&self, p: &QV) -> Vec<(usize, At)> {
        let x = qv_f64(p);
        let mut out = Vec::new();
        for (k, t) in self.tris.iter().enumerate() {
            if !t.near(&x) || t.plane(p).sign() != Ordering::Equal {
                continue;
            }
            if let Some(at) = t.at(p) {
                out.push((k, at));
            }
        }
        out
    }

    /// The crossings of the ray `p + s r` (`s > 0`) with the triangles not
    /// in `skip`, or `None` where the ray meets one's boundary or lies in its
    /// plane through it.
    fn crossings(&self, p: &QV, r: &[i64; 3], skip: &[usize]) -> Option<usize> {
        let x = qv_f64(p);
        let rf = r.map(|v| v as f64);
        let rv: V = r.map(ri);
        let mut count = 0;
        for (k, t) in self.tris.iter().enumerate() {
            if skip.contains(&k) {
                continue;
            }
            // The ray's parameters within the box (every direction positive).
            let (mut t0, mut t1) = (0.0f64, f64::INFINITY);
            for i in 0..3 {
                t0 = t0.max((t.lo[i] - x[i]) / rf[i]);
                t1 = t1.min((t.hi[i] - x[i]) / rf[i]);
            }
            if t0 > t1 + 1e-9 * (1.0 + t1.abs()) {
                continue;
            }
            let denom = sign(&dot(&t.n, &rv));
            // The plane's value at the point: the crossing ahead where it
            // and the direction's differ in sign.
            let s = t.plane(p).sign();
            if denom == Ordering::Equal {
                if s == Ordering::Equal {
                    return None;
                }
                continue;
            }
            // In the plane outside the triangle (`holding` took those inside
            // it), or the plane behind: no crossing ahead.
            if s == Ordering::Equal || s == denom {
                continue;
            }
            // The crossing point `p + s r` within each edge: `(n . r) (w .
            // x) - (w . r) (n . x)` over the point's offset, times the sign
            // of `n . r`.
            let nr = dot(&t.n, &rv);
            let mut inside = true;
            for k2 in 0..3 {
                let w = &t.w[k2];
                let g = sub(&scale(w, &nr), &scale(&t.n, &dot(w, &rv)));
                let v = qdot(&qsub(p, &qv(&t.c[k2])), &g).sign();
                match times(v, denom) {
                    Ordering::Equal => return None,
                    Ordering::Less => inside = false,
                    Ordering::Greater => {}
                }
            }
            if inside {
                count += 1;
            }
        }
        Some(count)
    }

    /// Where a point lies in the solid, pushed along directions in turn
    /// (symbolically: `p + e d1 + e^2 d2`, `e` infinitesimal).
    pub(super) fn member(&self, p: &QV, dirs: &[QV]) -> Loc {
        let held = self.holding(p);
        if held.is_empty() {
            for r in &RAYS {
                if let Some(n) = self.crossings(p, r, &[]) {
                    return if n % 2 == 1 { Loc::In } else { Loc::Out };
                }
            }
            return Loc::On;
        }
        if dirs.is_empty() {
            return Loc::On;
        }
        // On a triangle's wedge after the pushes: on the surface.
        for &(k, at) in &held {
            let t = &self.tris[k];
            if lex(dirs, &t.n) == Ordering::Equal
                && t.bounds(at).iter().all(|w| lex(dirs, w) != Ordering::Less)
            {
                return Loc::On;
            }
        }
        let skip: Vec<usize> = held.iter().map(|h| h.0).collect();
        'rays: for r in &RAYS {
            let rv: V = r.map(ri);
            let Some(far) = self.crossings(p, r, &skip) else {
                continue;
            };
            let mut near = 0;
            for &(k, at) in &held {
                let t = &self.tris[k];
                let nr = dot(&t.n, &rv);
                let snr = sign(&nr);
                if snr == Ordering::Equal {
                    continue 'rays;
                }
                // The pushed point's crossing with the wedge's plane lies
                // ahead where its side and the direction's differ.
                let nv = lex(dirs, &t.n);
                if nv == Ordering::Equal || nv == snr {
                    continue;
                }
                let mut inside = true;
                for w in t.bounds(at) {
                    let g = sub(&scale(w, &nr), &scale(&t.n, &dot(w, &rv)));
                    match times(lex(dirs, &g), snr) {
                        Ordering::Equal => continue 'rays,
                        Ordering::Less => inside = false,
                        Ordering::Greater => {}
                    }
                }
                if inside {
                    near += 1;
                }
            }
            return if (far + near) % 2 == 1 {
                Loc::In
            } else {
                Loc::Out
            };
        }
        Loc::On
    }

    /// Where a point on a face's (triangle's) plane lies in it.
    pub(super) fn in_face(&self, fi: usize, p: &QV) -> Loc {
        let t = &self.tris[fi];
        let mut on = false;
        for k in 0..3 {
            match t.side(k, p).sign() {
                Ordering::Less => return Loc::Out,
                Ordering::Equal => on = true,
                Ordering::Greater => {}
            }
        }
        if on {
            Loc::On
        } else {
            Loc::In
        }
    }
}

/// Whether a solid is an imported polyhedron (S9e.4b.2).
pub(super) fn is_mesh(s: &Solid) -> bool {
    crate::solid::boolean::polyhedra::imported::is_imported(s)
}

/// An imported polyhedron's model as an input `op`: each stored face's
/// triangles its faces on their exact planes, the stored edges and the
/// diagonals its edges, the stored vertices its vertices, named by the
/// stored ids.
pub(super) fn model(s: &Solid, op: Operand) -> Result<Prism> {
    let t = &s.topology;
    let triangles = crate::solid::boolean::polyhedra::imported::triangles(s)?;
    let id = |slot: Slot| {
        t.id_of(slot).ok_or(Error::InvalidTopology(
            "an imported polyhedron's unnamed slot",
        ))
    };
    let point = |p: crate::Point3| [p.x, p.y, p.z].map(crate::solid::split::q);
    // The stored vertices by their points, the stored edges by their ends.
    let mut vertex_of: BTreeMap<V, usize> = BTreeMap::new();
    for (vi, v) in t.vertices().iter().enumerate() {
        vertex_of.insert(point(v.position), vi);
    }
    let mut edge_of: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (ei, e) in t.edges().iter().enumerate() {
        let (Some(a), Some(b)) = (e.start, e.end) else {
            return Err(Error::InvalidTopology("an imported polyhedron's open edge"));
        };
        edge_of.insert((a.0.min(b.0), a.0.max(b.0)), ei);
    }
    let mut tris = Vec::new();
    let mut faces = Vec::new();
    let mut boxes = Vec::new();
    let mut corners: Vec<[usize; 3]> = Vec::new();
    for (fi, c) in triangles {
        let face = &t.faces()[fi];
        let outward = crate::solid::imported::outward(t, fi)
            .ok_or(Error::InvalidTopology("an imported polyhedron's face"))?;
        let crate::topology::Surface::Plane(frame) = &face.surface else {
            return Err(Error::InvalidTopology("an imported polyhedron's face"));
        };
        let sense = if frame.normal().dot(outward) > 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        let ks = [0, 1, 2].map(|k| vertex_of.get(&c[k]).copied());
        let [Some(a), Some(b), Some(cc)] = ks else {
            return Err(Error::InvalidTopology("an imported polyhedron's corner"));
        };
        corners.push([a, b, cc]);
        let tri = Tri::new(c);
        faces.push(MFace {
            kind: FaceKind::Facet(tris.len()),
            surf: Surf::Plane {
                p: tri.c[0].clone(),
                m: tri.n.clone(),
            },
            id: id(Slot::Face(FaceId(fi)))?,
            stored: face.surface.clone(),
            sense,
        });
        boxes.push((tri.lo, tri.hi));
        tris.push(tri);
    }
    // Each triangle's edges: left the triangle running it with its corners.
    let mut sides: BTreeMap<(usize, usize), [Option<usize>; 2]> = BTreeMap::new();
    for (k, c) in corners.iter().enumerate() {
        for j in 0..3 {
            let (a, b) = (c[j], c[(j + 1) % 3]);
            let key = (a.min(b), a.max(b));
            let slot = &mut sides.entry(key).or_insert([None, None])[usize::from(a > b)];
            if slot.is_some() {
                return Err(Error::InvalidTopology(
                    "an imported polyhedron's edge of several triangles",
                ));
            }
            *slot = Some(k);
        }
    }
    let mut verts: Vec<MVert> = Vec::new();
    let mut vert_of: BTreeMap<usize, usize> = BTreeMap::new();
    let mut edges: Vec<MEdge> = Vec::new();
    for (&(a, b), s) in &sides {
        let [Some(left), Some(right)] = *s else {
            return Err(Error::InvalidTopology("an imported polyhedron's open edge"));
        };
        let mut vertex = |v: usize, verts: &mut Vec<MVert>| -> Result<usize> {
            if let Some(&k) = vert_of.get(&v) {
                return Ok(k);
            }
            verts.push(MVert {
                p: qv(&point(t.vertices()[v].position)),
                id: Some(id(Slot::Vertex(VertexId(v)))?),
            });
            vert_of.insert(v, verts.len() - 1);
            Ok(verts.len() - 1)
        };
        let (va, vb) = (vertex(a, &mut verts)?, vertex(b, &mut verts)?);
        let (pa, pb) = (
            point(t.vertices()[a].position),
            point(t.vertices()[b].position),
        );
        let stored = edge_of
            .get(&(a, b))
            .map(|&e| id(Slot::Edge(EdgeId(e))))
            .transpose()?;
        edges.push(MEdge {
            kind: EdgeKind::Facet(edges.len()),
            curve: Crv::Line {
                p: qv(&pa),
                d: sub(&pb, &pa),
            },
            arc: None,
            start: va,
            end: vb,
            faces: [left, right],
            id: stored,
        });
    }
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    // Its one solid region (a cavity's void another region).
    let solid = t
        .regions()
        .iter()
        .position(|r| r.kind == crate::topology::RegionKind::Solid)
        .ok_or(Error::InvalidTopology("an imported polyhedron's region"))?;
    let region = id(Slot::Region(RegionId(solid)))?;
    info.insert(region, (op, Role::External));
    Ok(Prism {
        frame: Frame3::xy(),
        tolerance: s.resolution(),
        region,
        f: Affine::new(&Frame3::xy())?,
        lo: crate::solid::split::zero(),
        hi: crate::solid::split::zero(),
        bounds: Vec::new(),
        faces,
        edges,
        verts,
        info,
        boxes,
        ball: None,
        funnel: None,
        ring: None,
        given: None,
        hull: None,
        mesh: Some(Box::new(Mesh {
            tris,
            shells: t
                .regions()
                .iter()
                .filter(|r| r.kind == crate::topology::RegionKind::Solid)
                .map(|r| r.shells.len())
                .sum(),
        })),
    })
}
