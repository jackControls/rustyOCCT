//! S9c.1's exact prism models: the prism's affine frame on its stored axes
//! as rationals (`o + u x + v y + w n`), its profile as exact lines and
//! arcs (a full circle split into two arcs at a rational seam), its faces,
//! edges and vertices named by the input's slots, and exact tests of a
//! point against a face's region and against the solid, the latter with a
//! symbolic push along given directions (a point on the boundary decided
//! by where the push takes it).
use super::num::*;
use crate::identity::{EntityId, Role};
use crate::profile::boolean::Operand;
use crate::profile::{BoundaryKind, Segment};
use crate::solid::split::{q, zero};
use crate::solid::{Construction, Solid};
use crate::topology::{EdgeId, FaceId, Orientation, Slot, Surface, VertexId};
use crate::{Error, Frame3, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A frame's exact affine map on its stored axes.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Affine {
    pub(super) o: V,
    pub(super) x: V,
    pub(super) y: V,
    pub(super) n: V,
    /// The inverse's rows: local = inv . (p - o).
    inv: [V; 3],
}

impl Affine {
    pub(super) fn new(f: &Frame3) -> Result<Self> {
        let v = |a: [f64; 3]| a.map(q);
        let (o, x, y, n) = (
            v(f.origin().to_array()),
            v(f.x().to_array()),
            v(f.y().to_array()),
            v(f.normal().to_array()),
        );
        let det = dot(&x, &cross(&y, &n));
        if det == zero() {
            return Err(Error::Degenerate("a frame's axes"));
        }
        let inv = [
            scale(&cross(&y, &n), &(int(1) / &det)),
            scale(&cross(&n, &x), &(int(1) / &det)),
            scale(&cross(&x, &y), &(int(1) / &det)),
        ];
        Ok(Self { o, x, y, n, inv })
    }

    pub(super) fn point(&self, u: &R, v: &R, w: &R) -> V {
        add(
            &self.o,
            &add(
                &scale(&self.x, u),
                &add(&scale(&self.y, v), &scale(&self.n, w)),
            ),
        )
    }

    pub(super) fn vector(&self, u: &R, v: &R, w: &R) -> V {
        add(
            &scale(&self.x, u),
            &add(&scale(&self.y, v), &scale(&self.n, w)),
        )
    }

    pub(super) fn local(&self, p: &V) -> V {
        let d = sub(p, &self.o);
        [
            dot(&self.inv[0], &d),
            dot(&self.inv[1], &d),
            dot(&self.inv[2], &d),
        ]
    }

    pub(super) fn local_q(&self, p: &QV) -> QV {
        let d = qsub(p, &qv(&self.o));
        [
            qdot(&d, &self.inv[0]),
            qdot(&d, &self.inv[1]),
            qdot(&d, &self.inv[2]),
        ]
    }

    pub(super) fn local_dir(&self, d: &V) -> V {
        [
            dot(&self.inv[0], d),
            dot(&self.inv[1], d),
            dot(&self.inv[2], d),
        ]
    }

    pub(super) fn local_dir_q(&self, d: &QV) -> QV {
        [
            qdot(d, &self.inv[0]),
            qdot(d, &self.inv[1]),
            qdot(d, &self.inv[2]),
        ]
    }

    /// The covector of a local coordinate `k` in world terms (its gradient).
    pub(super) fn row(&self, k: usize) -> &V {
        &self.inv[k]
    }

    /// Whether the stored axes are exactly orthonormal: cylinders on them
    /// are circular.
    pub(super) fn orthonormal(&self) -> bool {
        let one = int(1);
        dot(&self.x, &self.x) == one
            && dot(&self.y, &self.y) == one
            && dot(&self.n, &self.n) == one
            && dot(&self.x, &self.y) == zero()
            && dot(&self.x, &self.n) == zero()
            && dot(&self.y, &self.n) == zero()
    }

    pub(super) fn same_axes(&self, o: &Self) -> bool {
        self.x == o.x && self.y == o.y && self.n == o.n
    }
}

pub(super) type P2 = [R; 2];

/// A profile segment, exact: its points are the stored binary64 values.
#[derive(Debug, Clone)]
pub(super) enum Seg {
    Line {
        p: P2,
        q: P2,
    },
    /// From `p` to `q` about `c` turning counter-clockwise (`ccw`) or not;
    /// both ends lie on the circle exactly.
    Arc {
        c: P2,
        r: R,
        p: P2,
        q: P2,
        ccw: bool,
    },
}

impl Seg {
    pub(super) fn start(&self) -> &P2 {
        match self {
            Seg::Line { p, .. } | Seg::Arc { p, .. } => p,
        }
    }
    pub(super) fn end(&self) -> &P2 {
        match self {
            Seg::Line { q, .. } | Seg::Arc { q, .. } => q,
        }
    }
    /// The unit-free tangent at a point of the segment, along its run.
    pub(super) fn tangent_at(&self, x: &[Qd; 2]) -> [Qd; 2] {
        match self {
            Seg::Line { p, q } => [Qd::rat(&q[0] - &p[0]), Qd::rat(&q[1] - &p[1])],
            Seg::Arc { c, ccw, .. } => {
                let (dx, dy) = (x[0].add_r(&-&c[0]), x[1].add_r(&-&c[1]));
                if *ccw {
                    [dy.neg(), dx]
                } else {
                    [dy, dx.neg()]
                }
            }
        }
    }
}

fn p2(p: crate::Point2) -> P2 {
    [q(p.x), q(p.y)]
}

/// A boundary of the profile: its segments and whether it is a hole
/// (material on its right: every boundary runs counter-clockwise).
#[derive(Debug, Clone)]
pub(super) struct Bound {
    pub(super) segs: Vec<Seg>,
    pub(super) hole: bool,
    /// A full circle split at a seam (its joints and verticals are not the
    /// input's).
    pub(super) circle: bool,
}

/// The rational point of a circle at the tangent of its half angle `s`.
pub(super) fn circle_point(c: &P2, r: &R, s: &R) -> P2 {
    let one = int(1);
    let den = &one + s * s;
    [
        &c[0] + r * (&one - s * s) / &den,
        &c[1] + r * (int(2) * s) / &den,
    ]
}

/// Where a result may cross the solid's boundary exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Loc {
    In,
    Out,
    /// On the boundary, even after every push.
    On,
}

/// A point of the profile against it: inside, outside, or on an element
/// (a segment's interior or a joint between two).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OnProfile {
    In,
    Out,
    Segment(usize, usize),
    Joint(usize, usize),
}

/// A face of a model: its kind and exact surface, and the input face it is
/// (its id and stored surface).
#[derive(Debug, Clone)]
pub(super) struct MFace {
    pub(super) kind: FaceKind,
    pub(super) surf: Surf,
    pub(super) id: EntityId,
    pub(super) stored: Surface,
    /// The stored face's sense: its surface's normal leaves the material
    /// when forward.
    pub(super) sense: Orientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum FaceKind {
    /// The cap at the high end or not.
    Cap(bool),
    /// The wall on segment `j` of boundary `b`.
    Wall(usize, usize),
    /// A sphere's hemisphere on the split plane's positive side or not
    /// (S9d.1).
    Half(bool),
    /// A cone's wall (S9d.3a).
    ConeWall,
    /// A torus's patch: the upper half of the tube (from its parallel
    /// seam) or not, the plus half of the turn (from its meridian seam) or
    /// not (S9d.4a).
    Patch(bool, bool),
}

/// A face's exact surface, its normal leaving the material.
#[derive(Debug, Clone)]
pub(super) enum Surf {
    Plane {
        p: V,
        m: V,
    },
    /// The prism's cylinder over a profile circle; the material inside the
    /// circle locally (`inside`) or outside it.
    Cyl {
        c: P2,
        r: R,
        inside: bool,
    },
    /// A sphere (S9d.1), its material inside.
    Sphere {
        c: V,
        r: R,
    },
    /// A cone's wall on the model's frame (S9d.3a): `u^2 + v^2 = (b + k
    /// w)^2`, its material inside.
    Cone {
        b: R,
        k: R,
    },
    /// A torus on the model's frame (S9d.4a), its material inside (its
    /// radii in the model's `ring`).
    Torus,
}

/// A 3D curve, exact.
#[derive(Debug, Clone)]
pub(super) enum Crv {
    /// `p + t d`.
    Line { p: QV, d: V },
    /// `c + a cos t + b sin t`.
    Conic { c: V, a: V, b: V },
    /// A piece of two cylinders' meeting (S9c.2), placed by its carrier's
    /// angle.
    Meet(Box<super::procedural::MeetCrv>),
    /// A circle of a surd radius (S9d.1: a sphere's sections and rims),
    /// placed by its basis coordinates.
    Circle(Box<super::sphere::Circ>),
    /// A cylinder's and a sphere's meeting over the cylinder's height
    /// (S9d.2b), placed by its height.
    Rise(Box<super::spheres::RiseCrv>),
    /// A plane's section of a cone (S9d.3a), placed by the cone's angle.
    Cone(Box<super::cone::ConeSec>),
    /// A plane's section of a torus (S9d.4a), a graph over one of its
    /// angles, placed by that angle's direction.
    Torus(Box<super::torus::TorusSec>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum EdgeKind {
    /// The cap edge (high end or not) over segment `j` of boundary `b`.
    Cap(bool, usize, usize),
    /// The vertical edge through the start of segment `j` of boundary `b`.
    Vertical(usize, usize),
    /// A sphere's rim (the high end or not), half `j` (S9d.1).
    Rim(bool, usize),
    /// A sphere's split great circle, arc `j` (S9d.1).
    Split(usize),
}

#[derive(Debug, Clone)]
pub(super) struct MEdge {
    pub(super) kind: EdgeKind,
    pub(super) curve: Crv,
    /// On a conic (or a circle), the start's and end's places and the
    /// direction.
    pub(super) arc: Option<([Qd; 2], [Qd; 2], bool)>,
    pub(super) start: usize,
    pub(super) end: usize,
    /// The faces on its left and right as it runs (model indices).
    pub(super) faces: [usize; 2],
    /// The input edge, none for a seam's.
    pub(super) id: Option<EntityId>,
}

#[derive(Debug, Clone)]
pub(super) struct MVert {
    pub(super) p: QV,
    pub(super) id: Option<EntityId>,
}

/// A prism's exact model.
#[derive(Debug, Clone)]
pub(super) struct Prism {
    /// The stored frame (arc frames and caps are built on it).
    pub(super) frame: Frame3,
    pub(super) tolerance: crate::Tolerance,
    /// The input's solid region.
    pub(super) region: EntityId,
    pub(super) f: Affine,
    pub(super) lo: R,
    pub(super) hi: R,
    pub(super) bounds: Vec<Bound>,
    pub(super) faces: Vec<MFace>,
    pub(super) edges: Vec<MEdge>,
    pub(super) verts: Vec<MVert>,
    /// Each input entity's operand and role.
    pub(super) info: BTreeMap<EntityId, (Operand, Role)>,
    /// Float bounds of each face, widened (for filtering only).
    pub(super) boxes: Vec<([f64; 3], [f64; 3])>,
    /// A sphere's own data (S9d.1): its faces are its ends' discs and its
    /// hemispheres.
    pub(super) ball: Option<super::sphere::Ball>,
    /// A cone's own data (S9d.3a): its faces are its ends' discs and its
    /// wall, over heights `lo..hi`.
    pub(super) funnel: Option<super::cone::Funnel>,
    /// A torus's own data (S9d.4a): its faces are its wall's four patches.
    pub(super) ring: Option<super::torus::Ring>,
}

fn out_of_domain(what: &'static str) -> Error {
    Error::OutOfDomain(what)
}

impl Prism {
    /// The exact model of a prism of lines and arcs; `seam` picks the
    /// rational point where full circles are split.
    pub(super) fn new(solid: &Solid, op: Operand, seam: &R) -> Result<Self> {
        let Construction::Prism(profile) = &solid.construction else {
            return Err(out_of_domain(
                "a Boolean of a solid other than a prism with arcs in any position (S9c)",
            ));
        };
        let f = Affine::new(&solid.frame)?;
        let (lo, hi) = (q(solid.start.min(solid.end)), q(solid.start.max(solid.end)));
        let t = &solid.topology;
        let id = |slot: Slot| {
            t.id_of(slot)
                .ok_or(Error::InvalidTopology("an unnamed slot"))
        };
        let mut bounds = Vec::new();
        for (b, boundary) in profile.boundaries().enumerate() {
            let hole = b > 0;
            let (segs, circle) = match &boundary.kind {
                BoundaryKind::Polygon(points) => (
                    (0..points.len())
                        .map(|j| Seg::Line {
                            p: p2(points[j]),
                            q: p2(points[(j + 1) % points.len()]),
                        })
                        .collect(),
                    false,
                ),
                BoundaryKind::Circle { center, radius } => {
                    let (c, r) = (p2(*center), q(*radius));
                    let s0 = circle_point(&c, &r, seam);
                    let s1 = [int(2) * &c[0] - &s0[0], int(2) * &c[1] - &s0[1]];
                    (
                        vec![
                            Seg::Arc {
                                c: c.clone(),
                                r: r.clone(),
                                p: s0.clone(),
                                q: s1.clone(),
                                ccw: true,
                            },
                            Seg::Arc {
                                c,
                                r,
                                p: s1,
                                q: s0,
                                ccw: true,
                            },
                        ],
                        true,
                    )
                }
                BoundaryKind::Path { points, segments } => {
                    let mut segs = Vec::new();
                    for j in 0..points.len() {
                        let (p, qq) = (p2(points[j]), p2(points[(j + 1) % points.len()]));
                        segs.push(match &segments[j] {
                            Segment::Line => Seg::Line { p, q: qq },
                            Segment::Arc {
                                center,
                                radius,
                                ccw,
                            } => {
                                let (c, r) = (p2(*center), q(*radius));
                                let on = |x: &P2| {
                                    let (dx, dy) = (&x[0] - &c[0], &x[1] - &c[1]);
                                    &dx * &dx + &dy * &dy == &r * &r
                                };
                                if !on(&p) || !on(&qq) {
                                    return Err(out_of_domain(
                                        "an arc whose ends lie off its circle in a Boolean of arcs in any position (S9c)",
                                    ));
                                }
                                Seg::Arc {
                                    c,
                                    r,
                                    p,
                                    q: qq,
                                    ccw: *ccw,
                                }
                            }
                            Segment::Spline(_) => {
                                return Err(out_of_domain(
                                    "a spline profile in a Boolean of prisms in any position (S9c)",
                                ))
                            }
                        });
                    }
                    (segs, false)
                }
            };
            bounds.push(Bound { segs, hole, circle });
        }
        // Faces: the caps, then the walls in profile order (the input's
        // slots: `Topology::prism`); a circle's halves share its wall.
        let n = f.n.clone();
        // A cap's exact plane holds the stored `x` and `y` (not normal to
        // the stored `n` unless the axes are orthonormal): its normal is
        // `x * y`, toward `n`.
        let up = {
            let xy = cross(&f.x, &f.y);
            if dot(&xy, &n) > zero() {
                xy
            } else {
                neg(&xy)
            }
        };
        let mut faces = vec![
            MFace {
                kind: FaceKind::Cap(false),
                surf: Surf::Plane {
                    p: f.point(&zero(), &zero(), &lo),
                    m: neg(&up),
                },
                id: id(Slot::Face(FaceId(0)))?,
                stored: t.faces()[0].surface.clone(),
                sense: t.faces()[0].sense,
            },
            MFace {
                kind: FaceKind::Cap(true),
                surf: Surf::Plane {
                    p: f.point(&zero(), &zero(), &hi),
                    m: up.clone(),
                },
                id: id(Slot::Face(FaceId(1)))?,
                stored: t.faces()[1].surface.clone(),
                sense: t.faces()[1].sense,
            },
        ];
        let mut wall_of: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        let mut next_face = 2;
        // Input vertex and edge slots, in `Topology::prism`'s order.
        let (mut vbase, mut ebase) = (0usize, 0usize);
        let mut verts: Vec<MVert> = Vec::new();
        let mut edges: Vec<MEdge> = Vec::new();
        for (b, bound) in bounds.iter().enumerate() {
            let count = bound.segs.len();
            let input_face = next_face;
            for (j, seg) in bound.segs.iter().enumerate() {
                let fid = if bound.circle {
                    input_face
                } else {
                    input_face + j
                };
                let surf = match seg {
                    Seg::Line { p, q: qq } => {
                        // Outward: the right of the run on an outer
                        // boundary, the left on a hole's.
                        let d = f.vector(&(&qq[0] - &p[0]), &(&qq[1] - &p[1]), &zero());
                        let m = cross(&d, &n);
                        Surf::Plane {
                            p: f.point(&p[0], &p[1], &lo),
                            m: if bound.hole { neg(&m) } else { m },
                        }
                    }
                    Seg::Arc { c, r, ccw, .. } => Surf::Cyl {
                        c: c.clone(),
                        r: r.clone(),
                        inside: *ccw != bound.hole,
                    },
                };
                wall_of.insert((b, j), faces.len());
                faces.push(MFace {
                    kind: FaceKind::Wall(b, j),
                    surf,
                    id: id(Slot::Face(FaceId(fid)))?,
                    stored: t.faces()[fid].surface.clone(),
                    sense: t.faces()[fid].sense,
                });
            }
            next_face += if bound.circle { 1 } else { count };
            // Vertices at both ends of each segment start.
            let first = verts.len();
            for (end, h) in [(false, &lo), (true, &hi)] {
                for (j, seg) in bound.segs.iter().enumerate() {
                    let p = seg.start();
                    let vid = if bound.circle {
                        None
                    } else {
                        let k = vbase + if end { count } else { 0 } + j;
                        Some(id(Slot::Vertex(VertexId(k)))?)
                    };
                    verts.push(MVert {
                        p: qv(&f.point(&p[0], &p[1], h)),
                        id: vid,
                    });
                }
            }
            let vat = |end: bool, j: usize| first + if end { count } else { 0 } + (j % count);
            let polygon = bound.segs.iter().all(|s| matches!(s, Seg::Line { .. }));
            for (j, seg) in bound.segs.iter().enumerate() {
                for (end, h) in [(false, &lo), (true, &hi)] {
                    let eid = if bound.circle {
                        Some(id(Slot::Edge(EdgeId(ebase + usize::from(end))))?)
                    } else if polygon {
                        Some(id(Slot::Edge(EdgeId(
                            ebase + if end { count } else { 0 } + j,
                        )))?)
                    } else {
                        Some(id(Slot::Edge(EdgeId(ebase + 2 * j + usize::from(end))))?)
                    };
                    let (curve, arc) = match seg {
                        Seg::Line { p, q: qq } => {
                            let a = f.point(&p[0], &p[1], h);
                            let b2 = f.point(&qq[0], &qq[1], h);
                            (
                                Crv::Line {
                                    p: qv(&a),
                                    d: sub(&b2, &a),
                                },
                                None,
                            )
                        }
                        Seg::Arc {
                            c,
                            r,
                            p,
                            q: qq,
                            ccw,
                        } => {
                            let rel = |x: &P2| [(&x[0] - &c[0]) / r, (&x[1] - &c[1]) / r];
                            (
                                Crv::Conic {
                                    c: f.point(&c[0], &c[1], h),
                                    a: scale(&f.x, r),
                                    b: scale(&f.y, r),
                                },
                                Some((rel(p).map(Qd::rat), rel(qq).map(Qd::rat), *ccw)),
                            )
                        }
                    };
                    // The cap's material lies on the run's left seen from
                    // outside the top cap (on an outer boundary): the wall
                    // is on its right there, the cap on its left.
                    let wall = wall_of[&(b, j)];
                    let cap = usize::from(end);
                    let left_is_cap = end != bound.hole;
                    edges.push(MEdge {
                        kind: EdgeKind::Cap(end, b, j),
                        curve,
                        arc,
                        start: vat(end, j),
                        end: vat(end, j + 1),
                        faces: if left_is_cap {
                            [cap, wall]
                        } else {
                            [wall, cap]
                        },
                        id: eid,
                    });
                }
            }
            for j in 0..count {
                let eid = if bound.circle {
                    None
                } else {
                    Some(id(Slot::Edge(EdgeId(ebase + 2 * count + j)))?)
                };
                let prev = wall_of[&(b, (j + count - 1) % count)];
                let this = wall_of[&(b, j)];
                let rational = |v: &MVert| {
                    v.p.clone()
                        .map(|x| x.rational().expect("a rational vertex").clone())
                };
                let a = rational(&verts[vat(false, j)]);
                let b2 = rational(&verts[vat(true, j)]);
                // Upward, the previous wall on the left (outside, facing
                // the wall's normal) on an outer boundary.
                edges.push(MEdge {
                    kind: EdgeKind::Vertical(b, j),
                    curve: Crv::Line {
                        p: qv(&a),
                        d: sub(&b2, &a),
                    },
                    arc: None,
                    start: vat(false, j),
                    end: vat(true, j),
                    faces: if bound.hole {
                        [this, prev]
                    } else {
                        [prev, this]
                    },
                    id: eid,
                });
            }
            if bound.circle {
                ebase += 2;
            } else {
                vbase += 2 * count;
                ebase += 3 * count;
            }
        }
        let mut info = BTreeMap::new();
        for (eid, _) in t.ids() {
            let role = t.derivation(eid).map_or(Role::External, |d| d.role);
            info.insert(eid, (op, role));
        }
        let mut prism = Self {
            frame: solid.frame,
            tolerance: solid.resolution(),
            region: id(Slot::Region(crate::topology::RegionId(1)))?,
            f,
            lo,
            hi,
            bounds,
            faces,
            edges,
            verts,
            info,
            boxes: Vec::new(),
            ball: None,
            funnel: None,
            ring: None,
        };
        prism.check_slots(solid)?;
        prism.boxes = (0..prism.faces.len()).map(|i| prism.face_box(i)).collect();
        Ok(prism)
    }

    /// The model's vertices against the input's stored ones (binary64): the
    /// slots mapped by `Topology::prism`'s order hold what they should.
    fn check_slots(&self, solid: &Solid) -> Result<()> {
        let t = &solid.topology;
        let tol = solid.resolution().linear();
        for v in &self.verts {
            let Some(vid) = v.id else { continue };
            let slot = t
                .ids()
                .find(|(i, _)| *i == vid)
                .map(|(_, s)| s)
                .ok_or(Error::InvalidTopology("a prism's vertex slot"))?;
            let Slot::Vertex(k) = slot else {
                return Err(Error::InvalidTopology("a prism's vertex slot"));
            };
            let stored = t.vertices()[k.0].position.to_array();
            let exact = qv_f64(&v.p);
            if (0..3).any(|i| (stored[i] - exact[i]).abs() > tol) {
                return Err(Error::InvalidTopology("a prism's vertices out of order"));
            }
        }
        Ok(())
    }

    /// A face's bounds in binary64, widened past rounding (for filtering).
    fn face_box(&self, fi: usize) -> ([f64; 3], [f64; 3]) {
        let f64s = |p: &V| p.clone().map(|x| crate::solid::split::rational_f64(&x));
        let mut pts: Vec<[f64; 3]> = Vec::new();
        let heights = [&self.lo, &self.hi];
        let add_seg = |seg: &Seg, pts: &mut Vec<[f64; 3]>| match seg {
            Seg::Line { p, q: qq } => {
                for h in heights {
                    pts.push(f64s(&self.f.point(&p[0], &p[1], h)));
                    pts.push(f64s(&self.f.point(&qq[0], &qq[1], h)));
                }
            }
            Seg::Arc { c, r, .. } => {
                for h in heights {
                    // The corners of the circle's square in the frame: its
                    // image holds the circle's (the frame may be turned).
                    for (dx, dy) in [(1, 1), (-1, 1), (1, -1), (-1, -1)] {
                        let u = &c[0] + r * int(dx);
                        let v = &c[1] + r * int(dy);
                        pts.push(f64s(&self.f.point(&u, &v, h)));
                    }
                }
            }
        };
        match self.faces[fi].kind {
            FaceKind::Cap(_) => {
                for b in &self.bounds {
                    for s in &b.segs {
                        add_seg(s, &mut pts);
                    }
                }
            }
            FaceKind::Wall(b, j) => add_seg(&self.bounds[b].segs[j], &mut pts),
            FaceKind::Half(_) | FaceKind::ConeWall | FaceKind::Patch(..) => {
                unreachable!("a sphere's, a cone's and a torus's boxes are their own")
            }
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &pts {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        for k in 0..3 {
            let m = 1e-9 * (1.0 + lo[k].abs().max(hi[k].abs()));
            lo[k] -= m;
            hi[k] += m;
        }
        (lo, hi)
    }

    /// The outward normal of a face at a point on it.
    pub(super) fn normal_at(&self, fi: usize, p: &QV) -> QV {
        match &self.faces[fi].surf {
            Surf::Plane { m, .. } => qv(m),
            Surf::Sphere { c, .. } => qsub(p, &qv(c)),
            Surf::Torus => {
                let ring = self.ring.as_ref().expect("a torus");
                let g = ring.gradient(&self.f.local_q(p));
                [0, 1, 2].map(|j| {
                    g[0].scale(&self.f.row(0)[j])
                        .add(&g[1].scale(&self.f.row(1)[j]))
                        .add(&g[2].scale(&self.f.row(2)[j]))
                })
            }
            Surf::Cone { b, k } => {
                // The gradient of u^2 + v^2 - (b + k w)^2 (halved).
                let l = self.f.local_q(p);
                let rk = l[2].scale(k).add_r(b).scale(k);
                [0, 1, 2].map(|j| {
                    l[0].scale(&self.f.row(0)[j])
                        .add(&l[1].scale(&self.f.row(1)[j]))
                        .sub(&rk.scale(&self.f.row(2)[j]))
                })
            }
            Surf::Cyl { c, inside, .. } => {
                let l = self.f.local_q(p);
                let (du, dv) = (l[0].add_r(&-&c[0]), l[1].add_r(&-&c[1]));
                // The gradient of (u - cu)^2 + (v - cv)^2: rows 0 and 1.
                let g = [0, 1, 2].map(|k| {
                    du.scale(&self.f.row(0)[k])
                        .add(&dv.scale(&self.f.row(1)[k]))
                });
                if *inside {
                    g
                } else {
                    g.map(|x| x.neg())
                }
            }
        }
    }

    // ------------------------------------------------------------ profile

    /// A profile point (of any field) against the profile, exactly: the
    /// parity of the chords crossed by the ray to `+u`, flipped inside each
    /// arc's circular segment (between its chord and itself). A point on an
    /// arc's chord counts as on the arc's side of it (the two tests jump
    /// together there).
    fn on_profile(&self, x: &[Qd; 2]) -> OnProfile {
        let mut inside = false;
        for (b, bound) in self.bounds.iter().enumerate() {
            for (j, seg) in bound.segs.iter().enumerate() {
                let (p, qq) = (seg.start(), seg.end());
                let xp = [x[0].add_r(&-&p[0]), x[1].add_r(&-&p[1])];
                if xp[0].sign() == Ordering::Equal && xp[1].sign() == Ordering::Equal {
                    return OnProfile::Joint(b, j);
                }
                let chord = [&qq[0] - &p[0], &qq[1] - &p[1]];
                // (q - p) x (x - p): positive on the chord's left.
                let mut side = xp[1].scale(&chord[0]).sub(&xp[0].scale(&chord[1])).sign();
                let between = || {
                    let t = xp[0].scale(&chord[0]).add(&xp[1].scale(&chord[1]));
                    let len = &chord[0] * &chord[0] + &chord[1] * &chord[1];
                    t.sign() == Ordering::Greater && t.cmp(&Qd::rat(len)) == Ordering::Less
                };
                match seg {
                    Seg::Line { .. } => {
                        if side == Ordering::Equal && between() {
                            return OnProfile::Segment(b, j);
                        }
                    }
                    Seg::Arc { c, r, ccw, .. } => {
                        let (dx, dy) = (x[0].add_r(&-&c[0]), x[1].add_r(&-&c[1]));
                        let g = dx.mul(&dx).add(&dy.mul(&dy)).add_r(&-(r * r)).sign();
                        if g == Ordering::Equal && self.within_arc(seg, &[dx, dy]) {
                            return OnProfile::Segment(b, j);
                        }
                        // A counter-clockwise arc lies on its chord's right.
                        let arc_side = if *ccw {
                            Ordering::Less
                        } else {
                            Ordering::Greater
                        };
                        if side == Ordering::Equal && between() {
                            side = arc_side;
                        }
                        if g == Ordering::Less && side == arc_side {
                            inside = !inside;
                        }
                    }
                }
                let below = |a: &R| x[1].cmp(&Qd::rat(a.clone())) == Ordering::Less;
                if below(&p[1]) != below(&qq[1]) {
                    // The ray crosses the chord where the point is on the
                    // chord's left for an upward chord.
                    let up = if qq[1] > p[1] {
                        Ordering::Greater
                    } else {
                        Ordering::Less
                    };
                    if side == up {
                        inside = !inside;
                    }
                }
            }
        }
        if inside {
            OnProfile::In
        } else {
            OnProfile::Out
        }
    }

    /// Whether a direction from an arc's centre lies strictly within its
    /// sweep.
    pub(super) fn within_arc(&self, seg: &Seg, d: &[Qd; 2]) -> bool {
        let Seg::Arc {
            c, p, q: qq, ccw, ..
        } = seg
        else {
            return false;
        };
        let pp = [&p[0] - &c[0], &p[1] - &c[1]];
        let qp = [&qq[0] - &c[0], &qq[1] - &c[1]];
        within_sweep(&pp, &qp, *ccw, d)
    }

    /// Where a profile point (of any field) lies, pushed along directions.
    fn profile_loc(&self, x: &[Qd; 2], dirs: &[[Qd; 2]]) -> Loc {
        match self.on_profile(x) {
            OnProfile::In => Loc::In,
            OnProfile::Out => Loc::Out,
            OnProfile::Segment(b, j) => match self.element_side(b, j, x, dirs) {
                Some(true) => Loc::In,
                Some(false) => Loc::Out,
                None => Loc::On,
            },
            OnProfile::Joint(b, j) => {
                let bound = &self.bounds[b];
                let count = bound.segs.len();
                let i = (j + count - 1) % count;
                let (inc, out) = (&bound.segs[i], &bound.segs[j]);
                let (ti, to) = (inc.tangent_at(x), out.tangent_at(x));
                // The turn at the joint: convex where the material is on the
                // turn's inside.
                let turn = ti[0].mul(&to[1]).sub(&ti[1].mul(&to[0])).sign();
                let a = self.element_side(b, i, x, dirs);
                let c = self.element_side(b, j, x, dirs);
                let left = !bound.hole;
                let r = match turn {
                    Ordering::Equal => {
                        // Smooth (or a cusp): the side of the element the push
                        // enters.
                        let Some(d) = dirs.first() else {
                            return Loc::On;
                        };
                        let along = to[0].mul(&d[0]).add(&to[1].mul(&d[1])).sign();
                        match along {
                            Ordering::Greater => c,
                            Ordering::Less => a,
                            Ordering::Equal => {
                                if a == c {
                                    a
                                } else {
                                    None
                                }
                            }
                        }
                    }
                    t => {
                        let convex = (t == Ordering::Greater) == left;
                        match (a, c) {
                            (Some(a), Some(c)) => Some(if convex { a && c } else { a || c }),
                            (Some(false), None) | (None, Some(false)) if convex => Some(false),
                            (Some(true), None) | (None, Some(true)) if !convex => Some(true),
                            _ => None,
                        }
                    }
                };
                match r {
                    Some(true) => Loc::In,
                    Some(false) => Loc::Out,
                    None => Loc::On,
                }
            }
        }
    }

    /// The material side of element `(b, j)`'s own curve at a point on it,
    /// pushed along directions (`None` when every push runs along it).
    fn element_side(&self, b: usize, j: usize, x: &[Qd; 2], dirs: &[[Qd; 2]]) -> Option<bool> {
        let bound = &self.bounds[b];
        match &bound.segs[j] {
            Seg::Line { p, q: qq } => {
                // Left of the run is material on an outer boundary.
                let d = [&qq[0] - &p[0], &qq[1] - &p[1]];
                let left = [-&d[1], d[0].clone()];
                for dir in dirs {
                    let s = dir[0].scale(&left[0]).add(&dir[1].scale(&left[1])).sign();
                    if s != Ordering::Equal {
                        return Some((s == Ordering::Greater) != bound.hole);
                    }
                }
                None
            }
            Seg::Arc { c, ccw, .. } => {
                // The side of the circle each push takes the point to, in
                // turn: the first-order change of |x - c|^2. A push along
                // the circle's tangent is one along the cylinder (a face on
                // it, or a face tangent to it, which is refused), not off it
                // by its curvature.
                let rel = [x[0].add_r(&-&c[0]), x[1].add_r(&-&c[1])];
                let orders: Vec<Qd> = dirs
                    .iter()
                    .map(|d| rel[0].mul(&d[0]).add(&rel[1].mul(&d[1])))
                    .collect();
                let material_inside = *ccw != bound.hole;
                for o in orders {
                    match o.sign() {
                        Ordering::Equal => continue,
                        s => return Some((s == Ordering::Less) == material_inside),
                    }
                }
                None
            }
        }
    }

    // ------------------------------------------------------------- solid

    /// Where a point lies in the solid, pushed along directions in turn
    /// (symbolically: `p + e d1 + e^2 d2`, `e` infinitesimal; a push along a
    /// cylinder's circle keeps to the cylinder).
    pub(super) fn member(&self, p: &QV, dirs: &[QV]) -> Loc {
        if let Some(ball) = &self.ball {
            return ball.member(p, dirs);
        }
        if let Some(fun) = &self.funnel {
            return fun.member(&self.f, &self.hi, p, dirs);
        }
        if let Some(ring) = &self.ring {
            return ring.member(&self.f, p, dirs);
        }
        let l = self.f.local_q(p);
        let ld: Vec<QV> = dirs.iter().map(|d| self.f.local_dir_q(d)).collect();
        // Heights.
        let mut height_in = true;
        for (bound, above) in [(&self.lo, true), (&self.hi, false)] {
            let mut s = l[2].add_r(&-bound).sign();
            let mut k = 0;
            while s == Ordering::Equal && k < ld.len() {
                s = ld[k][2].sign();
                k += 1;
            }
            let ok = match s {
                Ordering::Equal => return self.height_on(&l, &ld),
                Ordering::Greater => above,
                Ordering::Less => !above,
            };
            height_in &= ok;
        }
        let x = [l[0].clone(), l[1].clone()];
        let xd: Vec<[Qd; 2]> = ld.iter().map(|d| [d[0].clone(), d[1].clone()]).collect();
        let prof = self.profile_loc(&x, &xd);
        match (height_in, prof) {
            (false, _) => Loc::Out,
            (true, l) => l,
        }
    }

    /// A point on a cap's plane after every push: outside when off the
    /// profile, on the boundary otherwise.
    fn height_on(&self, l: &QV, ld: &[QV]) -> Loc {
        let x = [l[0].clone(), l[1].clone()];
        let xd: Vec<[Qd; 2]> = ld.iter().map(|d| [d[0].clone(), d[1].clone()]).collect();
        match self.profile_loc(&x, &xd) {
            Loc::Out => Loc::Out,
            _ => Loc::On,
        }
    }

    // ------------------------------------------------------------- faces

    /// Where a point on a face's surface lies in the face's region: inside,
    /// outside or on its boundary (exactly).
    pub(super) fn in_face(&self, fi: usize, p: &QV) -> Loc {
        if let Some(ball) = &self.ball {
            return match self.faces[fi].kind {
                FaceKind::Half(side) => ball.in_half(side, p),
                FaceKind::Cap(high) => {
                    let rim = ball.rim(usize::from(high)).expect("a disc's rim");
                    let d = qsub(p, &qv(&rim.c));
                    match qqdot(&d, &d).add_r(&-rim.r2.clone()).sign() {
                        Ordering::Less => Loc::In,
                        Ordering::Equal => Loc::On,
                        Ordering::Greater => Loc::Out,
                    }
                }
                FaceKind::Wall(..) | FaceKind::ConeWall | FaceKind::Patch(..) => {
                    unreachable!("a sphere has no walls")
                }
            };
        }
        if let Some(fun) = &self.funnel {
            return fun.in_face(&self.f, &self.hi, self.faces[fi].kind, p);
        }
        if let Some(ring) = &self.ring {
            return ring.in_face(&self.f, self.faces[fi].kind, p);
        }
        let l = self.f.local_q(p);
        let x = [l[0].clone(), l[1].clone()];
        match self.faces[fi].kind {
            FaceKind::Cap(_) => match self.on_profile(&x) {
                OnProfile::In => Loc::In,
                OnProfile::Out => Loc::Out,
                _ => Loc::On,
            },
            FaceKind::Half(_) | FaceKind::ConeWall | FaceKind::Patch(..) => {
                unreachable!("a prism has no hemispheres, cone walls or torus patches")
            }
            FaceKind::Wall(b, j) => {
                let h = [l[2].add_r(&-&self.lo).sign(), l[2].add_r(&-&self.hi).sign()];
                if h[0] == Ordering::Less || h[1] == Ordering::Greater {
                    return Loc::Out;
                }
                let on_height = h[0] == Ordering::Equal || h[1] == Ordering::Equal;
                let seg = &self.bounds[b].segs[j];
                let along = match seg {
                    Seg::Line { p, q: qq } => {
                        let d = [&qq[0] - &p[0], &qq[1] - &p[1]];
                        let t = x[0]
                            .add_r(&-&p[0])
                            .scale(&d[0])
                            .add(&x[1].add_r(&-&p[1]).scale(&d[1]));
                        let len = &d[0] * &d[0] + &d[1] * &d[1];
                        match (t.sign(), t.cmp(&Qd::rat(len))) {
                            (Ordering::Less, _) | (_, Ordering::Greater) => return Loc::Out,
                            (Ordering::Equal, _) | (_, Ordering::Equal) => Loc::On,
                            _ => Loc::In,
                        }
                    }
                    Seg::Arc { c, p, q: qq, .. } => {
                        let d = [x[0].add_r(&-&c[0]), x[1].add_r(&-&c[1])];
                        let at = |e: &P2| {
                            d[0].cmp(&Qd::rat(&e[0] - &c[0])) == Ordering::Equal
                                && d[1].cmp(&Qd::rat(&e[1] - &c[1])) == Ordering::Equal
                        };
                        // Directions only matter on the circle: compare the
                        // point's direction with the ends'.
                        if self.within_arc(seg, &d) {
                            Loc::In
                        } else if at(p) || at(qq) || on_ray(c, p, &d) || on_ray(c, qq, &d) {
                            Loc::On
                        } else {
                            return Loc::Out;
                        }
                    }
                };
                if on_height {
                    Loc::On
                } else {
                    along
                }
            }
        }
    }
}

/// Whether the direction `d` from `c` runs along the ray to `e`.
fn on_ray(c: &P2, e: &P2, d: &[Qd; 2]) -> bool {
    let r = [&e[0] - &c[0], &e[1] - &c[1]];
    let cr = d[0].scale(&r[1]).sub(&d[1].scale(&r[0])).sign();
    let dt = d[0].scale(&r[0]).add(&d[1].scale(&r[1])).sign();
    cr == Ordering::Equal && dt == Ordering::Greater
}

/// The pseudo-angle order of two directions, counter-clockwise from +x:
/// `Less` when `a` comes first.
pub(super) fn angle_cmp(a: &[Qd; 2], b: &[Qd; 2]) -> Ordering {
    let half = |v: &[Qd; 2]| {
        let s = v[1].sign();
        s == Ordering::Less || (s == Ordering::Equal && v[0].sign() == Ordering::Less)
    };
    match (half(a), half(b)) {
        (false, true) => Ordering::Less,
        (true, false) => Ordering::Greater,
        _ => {
            // Same half: b to the left of a comes later.
            let c = mixed_dot_sign(&[a[0].clone(), a[1].neg()], &[b[1].clone(), b[0].clone()]);
            match c {
                Ordering::Greater => Ordering::Less,
                Ordering::Less => Ordering::Greater,
                Ordering::Equal => Ordering::Equal,
            }
        }
    }
}

/// A direction relative to a rational base: its coordinates in the base's
/// rotation (`(d . e, e x d)`), mirrored for a clockwise sweep.
pub(super) fn relative(e: &P2, d: &[Qd; 2], ccw: bool) -> [Qd; 2] {
    let x = d[0].scale(&e[0]).add(&d[1].scale(&e[1]));
    let y = d[1].scale(&e[0]).sub(&d[0].scale(&e[1]));
    [x, if ccw { y } else { y.neg() }]
}

/// Whether `d` lies strictly within the sweep from `p` to `q` (directions
/// from the centre) turning counter-clockwise or not.
pub(super) fn within_sweep(p: &P2, q: &P2, ccw: bool, d: &[Qd; 2]) -> bool {
    let rd = relative(p, d, ccw);
    let rq = relative(p, &[Qd::rat(q[0].clone()), Qd::rat(q[1].clone())], ccw);
    let zero_dir = [Qd::rat(int(1)), Qd::rat(zero())];
    // The end at a full turn when q = p.
    let full = rq[1].sign() == Ordering::Equal && rq[0].sign() == Ordering::Greater;
    angle_cmp(&zero_dir, &rd) == Ordering::Less && (full || angle_cmp(&rd, &rq) == Ordering::Less)
}
