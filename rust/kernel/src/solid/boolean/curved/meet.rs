//! S9c.1's meetings, exact: an edge's points on the other solid's face
//! surfaces (a line or an arc against a plane or a cylinder), and the
//! curves where two faces' surfaces meet (lines, generatrices, plane
//! sections of cylinders, the two ellipses of equal cylinders whose axes
//! cross). Two cylinders meet in S9c.1's curves only when both are circular
//! (their frames' axes equal or exactly orthonormal); others are S9c.2's
//! unless their surfaces are certainly apart.
use super::model::{Affine, Crv, Prism, Surf, P2};
use super::num::*;
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// Where an edge meets a surface.
#[derive(Debug, Clone)]
pub(super) enum EdgeMeet {
    None,
    /// Points, each with its place along the edge's curve.
    Points(Vec<(Pos, QV)>),
    /// The edge's curve lies on the surface.
    Along,
}

/// A place along a curve: a line's parameter, or a conic's `(cos, sin)`.
#[derive(Debug, Clone)]
pub(super) enum Pos {
    T(Qd),
    Ang([Qd; 2]),
}

pub(super) fn tangency() -> Error {
    Error::Degenerate("a tangency between the inputs (S9c)")
}

fn quartic() -> Error {
    Error::OutOfDomain("two cylinders meeting in curves other than lines and conics (S9c.2)")
}

/// Solves `a t^2 + b t + c = 0` (`a != 0`): no root, or two (a double
/// root is a tangency).
pub(super) fn quadratic(a: &R, b: &R, c: &R) -> Result<Vec<Qd>> {
    let disc = b * b - int(4) * a * c;
    match sign(&disc) {
        Ordering::Less => Ok(Vec::new()),
        Ordering::Equal => Err(tangency()),
        Ordering::Greater => {
            let den = int(2) * a;
            let (x, y) = (-b / &den, int(1) / &den);
            Ok(vec![
                Qd::new(x.clone(), -&y, disc.clone()),
                Qd::new(x, y, disc),
            ])
        }
    }
}

/// Solves `alpha cos + beta sin = gamma` for `(cos, sin)`: none, or two.
pub(super) fn trig(alpha: &R, beta: &R, gamma: &R) -> Result<Option<Vec<[Qd; 2]>>> {
    let rho2 = alpha * alpha + beta * beta;
    if rho2 == zero() {
        // Every angle or none.
        return Ok(if *gamma == zero() {
            None
        } else {
            Some(Vec::new())
        });
    }
    let d = &rho2 - gamma * gamma;
    match sign(&d) {
        Ordering::Less => Ok(Some(Vec::new())),
        Ordering::Equal => Err(tangency()),
        Ordering::Greater => {
            let inv = int(1) / &rho2;
            let mut out = Vec::new();
            for s in [-1, 1] {
                let k = int(s);
                out.push([
                    Qd::new(alpha * gamma * &inv, -(beta * &k) * &inv, d.clone()),
                    Qd::new(beta * gamma * &inv, alpha * &k * &inv, d.clone()),
                ]);
            }
            Ok(Some(out))
        }
    }
}

/// A point of a conic at `(cos, sin)`.
pub(super) fn conic_point(c: &V, a: &V, b: &V, cs: &[Qd; 2]) -> QV {
    qadd(&qadd(&qv(c), &qscale(a, &cs[0])), &qscale(b, &cs[1]))
}

/// The tangent of a conic at `(cos, sin)` (increasing angle).
pub(super) fn conic_tangent(a: &V, b: &V, cs: &[Qd; 2]) -> QV {
    qadd(&qscale(a, &cs[1].neg()), &qscale(b, &cs[0]))
}

/// A plane's coefficients on a prism's cylinder: `alpha cos + beta sin +
/// mu w + kappa = 0` over `(cu + r cos, cv + r sin, w)`.
fn plane_on_cylinder(f: &Affine, c: &P2, r: &R, p: &V, m: &V) -> [R; 4] {
    let (mx, my, mn) = (dot(m, &f.x), dot(m, &f.y), dot(m, &f.n));
    let kappa = dot(m, &sub(&f.o, p)) + &c[0] * &mx + &c[1] * &my;
    [r * &mx, r * &my, mn, kappa]
}

/// A curve where two faces' surfaces meet: its branches (each a line or a
/// conic), or the surfaces coincide.
#[derive(Debug, Clone)]
pub(super) enum Section {
    Curves(Vec<Crv>),
    Same,
}

/// Equal circular cylinders with crossing axes: the two planes their
/// ellipses lie in, and the points where those cross.
#[derive(Debug, Clone)]
pub(super) struct Crossing {
    pub(super) planes: [(V, V); 2],
    pub(super) points: Vec<QV>,
}

/// How two cylinders (of different prisms) meet.
#[derive(Debug, Clone)]
pub(super) enum CylPair {
    Apart,
    /// Parallel circular cylinders: their circles in the first one's frame
    /// (`(u, v)` at any height).
    Parallel {
        c2: P2,
        r2: R,
    },
    /// Equal circular cylinders with crossing axes: the two planes their
    /// ellipses lie in, and the points where those cross.
    Crossing(Box<Crossing>),
    /// Circular cylinders with perpendicular axes meeting in a quartic
    /// (S9c.2a).
    Quartic(Box<super::procedural::Quartic>),
    /// Parallel cylinders not circular in a common measure (S9c.2b.2):
    /// generatrices along `d` through algebraic points.
    Lines(Vec<QV>, V),
    /// A cylinder and a sphere in loops (S9d.2b): pieces over the angle and
    /// the height, and their switches.
    Mixed(Box<super::spheres::Mixed>),
    /// Two cones whose quadrics differ by an affine function (S9d.3b):
    /// their meeting on the plane through `p` normal to `m`, the section of
    /// operand `k`'s cone.
    Plane {
        k: usize,
        p: V,
        m: V,
    },
    /// One surface.
    Same,
}

/// Two cylinders' relation: `x` a cylinder of `px`, `y` one of `py`.
pub(super) fn cyl_pair(
    px: &Prism,
    cx: &P2,
    rx: &R,
    boxes: [&([f64; 3], [f64; 3]); 2],
    py: &Prism,
    cy: &P2,
    ry: &R,
) -> Result<CylPair> {
    let (fx, fy) = (&px.f, &py.f);
    let apart_boxes =
        (0..3).any(|k| boxes[0].1[k] < boxes[1].0[k] || boxes[1].1[k] < boxes[0].0[k]);
    let circular = fx.same_axes(fy) || (fx.orthonormal() && fy.orthonormal());
    let parallel = is_zero(&cross(&fx.n, &fy.n));
    if parallel {
        // The other axis in the first frame: a point and whether circles.
        let oy = fy.point(&cy[0], &cy[1], &zero());
        let l = fx.local(&oy);
        if circular {
            // In the first frame the other circle has the same radius when
            // the axes are equal or both orthonormal (lengths kept).
            let c2 = [l[0].clone(), l[1].clone()];
            let d = [&c2[0] - &cx[0], &c2[1] - &cx[1]];
            let dist2 = &d[0] * &d[0] + &d[1] * &d[1];
            if dist2 == zero() {
                return Ok(if rx == ry {
                    CylPair::Same
                } else {
                    CylPair::Apart
                });
            }
            let (sum, diff) = (rx + ry, rx - ry);
            let (s2, d2) = (&sum * &sum, &diff * &diff);
            return match (dist2.cmp(&s2), dist2.cmp(&d2)) {
                (Ordering::Greater, _) | (_, Ordering::Less) => Ok(CylPair::Apart),
                (Ordering::Equal, _) | (_, Ordering::Equal) => Err(tangency()),
                _ => Ok(CylPair::Parallel { c2, r2: ry.clone() }),
            };
        }
        if apart_boxes || parallel_apart(fx, cx, rx, fy, cy, ry) {
            return Ok(CylPair::Apart);
        }
        return super::algebraic::parallel((fx, cx, rx), (fy, cy, ry));
    }
    if apart_boxes {
        return Ok(CylPair::Apart);
    }
    if !circular {
        // Turned frames, crossing axes (S9c.2b.1).
        return super::turned::crossing((fx, cx, rx), (fy, cy, ry), px.tolerance.linear());
    }
    // Circular, equal radii: the axes must meet (coplanar) for two conics;
    // others perpendicular (in exact frames) meet in quartics (S9c.2a).
    let ox = fx.point(&cx[0], &cx[1], &zero());
    let oy = fy.point(&cy[0], &cy[1], &zero());
    let w = cross(&fx.n, &fy.n);
    if rx != ry || dot(&sub(&oy, &ox), &w) != zero() {
        if dot(&fx.n, &fy.n) != zero() {
            return Err(quartic());
        }
        return super::procedural::perpendicular((fx, cx, rx), (fy, cy, ry));
    }
    // The axes' meeting point: ox + s nx = oy + t ny.
    let d = sub(&oy, &ox);
    let s = dot(&cross(&d, &fy.n), &w) / dot(&w, &w);
    let p = add(&ox, &scale(&fx.n, &s));
    // Unit axes (orthonormal frames): the bisector planes' normals.
    let planes = [
        (p.clone(), add(&fx.n, &fy.n)),
        (p.clone(), sub(&fx.n, &fy.n)),
    ];
    // Where the two ellipses cross: p +- r w / |w|.
    let w2 = dot(&w, &w);
    let k = Qd::new(zero(), rx / &w2, w2.clone());
    let points = [k.clone(), k.neg()]
        .iter()
        .map(|k| qadd(&qv(&p), &qscale(&w, k)))
        .collect();
    Ok(CylPair::Crossing(Box::new(Crossing { planes, points })))
}

/// Whether two parallel cylinders (not circular in a common frame) are
/// certainly apart: in the first frame's `(u, v)`, the other's section is
/// an ellipse whose points all lie within or beyond the first circle by a
/// certified margin.
fn parallel_apart(fx: &Affine, cx: &P2, rx: &R, fy: &Affine, cy: &P2, ry: &R) -> bool {
    use crate::certified::Interval as I;
    // The other's points in the first frame: l(cy + ry e) projected on (u, v).
    let oy = fy.point(&cy[0], &cy[1], &zero());
    let l = fx.local(&oy);
    let (ax, ay) = (
        fx.local_dir(&scale(&fy.x, ry)),
        fx.local_dir(&scale(&fy.y, ry)),
    );
    // Distance from the first centre: |l - cx| and the ellipse's reach
    // at most |a| + |b| in each coordinate.
    let d = [&l[0] - &cx[0], &l[1] - &cx[1]];
    let reach = |i: usize| &ax[i] * &ax[i] + &ay[i] * &ay[i];
    let norm = |x: &R| I::exact(x.clone());
    let dist = norm(&(&d[0] * &d[0] + &d[1] * &d[1])).sqrt();
    let spread = norm(&(reach(0) + reach(1))).sqrt();
    let r = I::exact(rx.clone());
    // Inside: dist + spread < r; outside: dist - spread > r.
    let inside = r.sub(&dist.add(&spread)).sign() == Some(Ordering::Greater);
    let outside = dist.sub(&spread).sub(&r).sign() == Some(Ordering::Greater);
    inside || outside
}

/// Where two faces' surfaces meet (face `fx` of `px`, `fy` of `py`).
pub(super) fn section(
    px: &Prism,
    fx: usize,
    py: &Prism,
    fy: usize,
    pair: Option<&CylPair>,
) -> Result<Section> {
    match (&px.faces[fx].surf, &py.faces[fy].surf) {
        // S9d.4a: a plane's section of a torus, by the pair's relation.
        (Surf::Plane { .. }, Surf::Torus) | (Surf::Torus, Surf::Plane { .. }) => {
            match pair.expect("a plane and a torus's relation") {
                CylPair::Mixed(x) => Ok(Section::Curves(x.pieces.clone())),
                _ => Ok(Section::Curves(Vec::new())),
            }
        }
        (Surf::Torus, _) | (_, Surf::Torus) => {
            Err(Error::OutOfDomain("a torus against a curved face (S9d.4b)"))
        }
        // S9d.3a: a plane's section of a cone's wall.
        (Surf::Plane { p, m }, Surf::Cone { .. }) => Ok(Section::Curves(super::cone::plane_cone(
            1,
            &py.f,
            py.funnel.as_ref().expect("a cone"),
            p,
            m,
        )?)),
        (Surf::Cone { .. }, Surf::Plane { p, m }) => Ok(Section::Curves(super::cone::plane_cone(
            0,
            &px.f,
            px.funnel.as_ref().expect("a cone"),
            p,
            m,
        )?)),
        // S9d.3b: a cone against a cylinder, a sphere or a cone, by their
        // relation (rings over a carrier, a plane, or apart).
        (Surf::Cone { .. }, _) | (_, Surf::Cone { .. }) => {
            match pair.expect("a cone and a curved face's relation") {
                CylPair::Quartic(x) => Ok(Section::Curves(
                    x.pieces
                        .iter()
                        .map(|m| Crv::Meet(Box::new(m.clone())))
                        .collect(),
                )),
                CylPair::Plane { k, p, m } => {
                    let cone = if *k == 0 { px } else { py };
                    Ok(Section::Curves(super::cone::plane_cone(
                        *k,
                        &cone.f,
                        cone.funnel.as_ref().expect("a cone"),
                        p,
                        m,
                    )?))
                }
                // Loops of a cone and a sphere (S9d.3b.2).
                CylPair::Mixed(x) => Ok(Section::Curves(x.pieces.clone())),
                CylPair::Same => Ok(Section::Same),
                CylPair::Apart => Ok(Section::Curves(Vec::new())),
                CylPair::Parallel { .. } | CylPair::Crossing(_) | CylPair::Lines(..) => {
                    unreachable!("two cylinders' relations")
                }
            }
        }
        (Surf::Plane { p: p1, m: m1 }, Surf::Plane { p: p2, m: m2 }) => {
            let d = cross(m1, m2);
            if is_zero(&d) {
                return Ok(if dot(m1, &sub(p2, p1)) == zero() {
                    Section::Same
                } else {
                    Section::Curves(Vec::new())
                });
            }
            let (k1, k2) = (dot(m1, p1), dot(m2, p2));
            let dd = dot(&d, &d);
            let p0 = scale(
                &add(&scale(&cross(m2, &d), &k1), &scale(&cross(&d, m1), &k2)),
                &(int(1) / dd),
            );
            Ok(Section::Curves(vec![Crv::Line { p: qv(&p0), d }]))
        }
        (Surf::Plane { p, m }, Surf::Cyl { c, r, .. }) => plane_cyl(&py.f, c, r, p, m),
        (Surf::Cyl { c, r, .. }, Surf::Plane { p, m }) => plane_cyl(&px.f, c, r, p, m),
        // A plane's section of a sphere: a circle (S9d.1).
        (Surf::Plane { p, m }, Surf::Sphere { c, r })
        | (Surf::Sphere { c, r }, Surf::Plane { p, m }) => Ok(Section::Curves(
            super::sphere::plane_section(c, r, p, m)?
                .map(|circ| Crv::Circle(Box::new(circ)))
                .into_iter()
                .collect(),
        )),
        // S9d.2: a cylinder and a sphere by their relation (rings over the
        // cylinder's angle, or apart); two spheres on their radical plane.
        (Surf::Cyl { .. }, Surf::Sphere { .. }) | (Surf::Sphere { .. }, Surf::Cyl { .. }) => {
            match pair.expect("a cylinder and a sphere's relation") {
                CylPair::Quartic(x) => Ok(Section::Curves(
                    x.pieces
                        .iter()
                        .map(|m| Crv::Meet(Box::new(m.clone())))
                        .collect(),
                )),
                CylPair::Mixed(x) => Ok(Section::Curves(x.pieces.clone())),
                _ => Ok(Section::Curves(Vec::new())),
            }
        }
        (Surf::Sphere { c: c1, r: r1 }, Surf::Sphere { c: c2, r: r2 }) => Ok(Section::Curves(
            super::spheres::sphere_sphere(c1, r1, c2, r2)?
                .map(|circ| Crv::Circle(Box::new(circ)))
                .into_iter()
                .collect(),
        )),
        (Surf::Cyl { c, r, .. }, Surf::Cyl { .. }) => {
            match pair.expect("a cylinder pair's relation") {
                CylPair::Mixed(_) | CylPair::Plane { .. } => {
                    unreachable!("a cylinder and a sphere's or two cones'")
                }
                CylPair::Apart => Ok(Section::Curves(Vec::new())),
                CylPair::Same => Ok(Section::Same),
                CylPair::Parallel { c2, r2 } => {
                    // The circles' meetings in the first frame: generatrices.
                    let cs = circle_circle(c, r, c2, r2)?;
                    let f = &px.f;
                    Ok(Section::Curves(
                        cs.iter()
                            .map(|e| {
                                let u = e[0].scale(r).add_r(&c[0]);
                                let v = e[1].scale(r).add_r(&c[1]);
                                let base =
                                    qadd(&qadd(&qv(&f.o), &qscale(&f.x, &u)), &qscale(&f.y, &v));
                                Crv::Line {
                                    p: base,
                                    d: f.n.clone(),
                                }
                            })
                            .collect(),
                    ))
                }
                CylPair::Lines(points, d) => Ok(Section::Curves(
                    points
                        .iter()
                        .map(|p| Crv::Line {
                            p: p.clone(),
                            d: d.clone(),
                        })
                        .collect(),
                )),
                CylPair::Quartic(x) => Ok(Section::Curves(
                    x.pieces
                        .iter()
                        .map(|m| Crv::Meet(Box::new(m.clone())))
                        .collect(),
                )),
                CylPair::Crossing(x) => {
                    let mut out = Vec::new();
                    for (p, m) in &x.planes {
                        if let Section::Curves(c2) = plane_cyl(&px.f, c, r, p, m)? {
                            out.extend(c2);
                        }
                    }
                    Ok(Section::Curves(out))
                }
            }
        }
    }
}

/// A plane's section of a prism's cylinder.
fn plane_cyl(f: &Affine, c: &P2, r: &R, p: &V, m: &V) -> Result<Section> {
    let [alpha, beta, mu, kappa] = plane_on_cylinder(f, c, r, p, m);
    if mu == zero() {
        // Parallel to the axis: generatrices (or none, or the surface
        // itself never: a plane is not a cylinder).
        let Some(roots) = trig(&alpha, &beta, &-&kappa)? else {
            return Ok(Section::Curves(Vec::new()));
        };
        return Ok(Section::Curves(
            roots
                .iter()
                .map(|e| {
                    let u = e[0].scale(r).add_r(&c[0]);
                    let v = e[1].scale(r).add_r(&c[1]);
                    Crv::Line {
                        p: qadd(&qadd(&qv(&f.o), &qscale(&f.x, &u)), &qscale(&f.y, &v)),
                        d: f.n.clone(),
                    }
                })
                .collect(),
        ));
    }
    // A plane within rounding of the axis's direction (its section an
    // ellipse whose axis along the cylinder is past 10^12 radii, beyond any
    // binary64 edge): the height's turn rate against the radius.
    if &alpha * &alpha + &beta * &beta > int(10).pow(24) * &mu * &mu * r * r {
        return Err(Error::Degenerate(
            "a plane within rounding of a cylinder's direction",
        ));
    }
    // w = -(kappa + alpha cos + beta sin) / mu.
    let inv = int(-1) / &mu;
    let c0 = f.point(&c[0], &c[1], &(&kappa * &inv));
    let a = add(&scale(&f.x, r), &scale(&f.n, &(&alpha * &inv)));
    let b = add(&scale(&f.y, r), &scale(&f.n, &(&beta * &inv)));
    Ok(Section::Curves(vec![Crv::Conic { c: c0, a, b }]))
}

/// Two circles' meetings in a plane: `(cos, sin)` on the first.
fn circle_circle(c1: &P2, r1: &R, c2: &P2, r2: &R) -> Result<Vec<[Qd; 2]>> {
    // |c1 + r1 e - c2|^2 = r2^2: 2 r1 (c1 - c2) . e = r2^2 - r1^2 - |c1 - c2|^2.
    let d = [&c1[0] - &c2[0], &c1[1] - &c2[1]];
    let alpha = int(2) * r1 * &d[0];
    let beta = int(2) * r1 * &d[1];
    let gamma = r2 * r2 - r1 * r1 - (&d[0] * &d[0] + &d[1] * &d[1]);
    Ok(trig(&alpha, &beta, &gamma)?.unwrap_or_default())
}

/// Where an edge's curve meets a face's surface (face `fy` of `py`).
pub(super) fn edge_surface(
    curve: &Crv,
    own: Option<(&Prism, &P2, &R)>,
    own_ball: Option<&super::sphere::Ball>,
    py: &Prism,
    fy: usize,
    pair: Option<&CylPair>,
) -> Result<EdgeMeet> {
    match (curve, &py.faces[fy].surf) {
        // S9d.4a: a line against a torus; a circle is S9d.4b's.
        (Crv::Line { p, d }, Surf::Torus) => {
            super::torus::line_torus(p, d, &py.f, py.ring.as_ref().expect("a torus"))
        }
        (Crv::Conic { .. } | Crv::Circle(_), Surf::Torus) => {
            Err(Error::OutOfDomain("a torus against a curved face (S9d.4b)"))
        }
        // S9d.3a: a line against a cone's wall; S9d.3b's others.
        (Crv::Line { p, d }, Surf::Cone { .. }) => {
            super::cone::line_cone(p, d, &py.f, py.funnel.as_ref().expect("a cone"), &py.hi)
        }
        // S9d.3b: a cap's or a rim's circle against a cone, a sphere's circle
        // against a cone.
        (Crv::Conic { c, a, b }, Surf::Cone { b: rb, k }) => {
            match super::algebraic::circle_quadric(
                c,
                a,
                b,
                &super::procedural::other_cone(&py.f, rb, k),
            ) {
                Ok(points) => Ok(EdgeMeet::Points(points)),
                // A circle on the cone's quadric past its ends: no meeting.
                Err(Error::Degenerate(_)) if beyond_ends(&py.f, &py.hi, c, a, b) => {
                    Ok(EdgeMeet::None)
                }
                Err(e) => Err(e),
            }
        }
        (Crv::Circle(circ), Surf::Cone { b: rb, k }) => {
            super::spheres::circ_quadric(circ, &super::procedural::other_cone(&py.f, rb, k))
        }
        // S9d.1: a line against a sphere; a sphere's circle against a
        // plane; anything else against a sphere is S9d.2's.
        (Crv::Line { p, d }, Surf::Sphere { c, r }) => super::sphere::line_sphere(p, d, c, r),
        (Crv::Circle(circ), Surf::Plane { p: p0, m }) => match circ.meet_plane(p0, m)? {
            None => Ok(EdgeMeet::Along),
            Some(cs) => Ok(EdgeMeet::Points(
                cs.into_iter()
                    .map(|e| {
                        let x = qadd(
                            &qv(&circ.c),
                            &qadd(&qscale(&circ.x, &e[0]), &qscale(&circ.y, &e[1])),
                        );
                        (Pos::Ang(e), x)
                    })
                    .collect(),
            )),
        },
        // S9d.2: a sphere's circle against a cylinder or another sphere, a
        // prism's arc against a sphere.
        (Crv::Circle(circ), Surf::Cyl { c: cy, r: ry, .. }) => {
            super::spheres::circ_cylinder(circ, &py.f, cy, ry)
        }
        (Crv::Circle(circ), Surf::Sphere { c: c2, r: r2 }) => {
            let ball = own_ball.expect("a sphere's circle's own sphere");
            super::spheres::circ_sphere(circ, &ball.c, &ball.r, c2, r2)
        }
        (Crv::Conic { c, a, b }, Surf::Sphere { c: cs, r }) => Ok(EdgeMeet::Points(
            super::algebraic::circle_quadric(c, a, b, &super::procedural::other_sphere(cs, r))?,
        )),
        (Crv::Meet(_) | Crv::Rise(_) | Crv::Cone(_) | Crv::Torus(_), _) => {
            unreachable!("a model edge is a line, an arc or a circle")
        }
        (Crv::Line { p, d }, Surf::Plane { p: p0, m }) => {
            let md = dot(m, d);
            let off = qdot(&qsub(&qv(p0), p), m);
            if md == zero() {
                return Ok(if off.sign() == Ordering::Equal {
                    EdgeMeet::Along
                } else {
                    EdgeMeet::None
                });
            }
            let t = off.scale(&(int(1) / md));
            let x = qadd(p, &qscale(d, &t));
            Ok(EdgeMeet::Points(vec![(Pos::T(t), x)]))
        }
        (Crv::Line { p, d }, Surf::Cyl { c, r, .. }) => {
            let f = &py.f;
            let l = f.local_q(p);
            let ld = f.local_dir(d);
            // (lu + t du - cu)^2 + (lv + t dv - cv)^2 = r^2.
            let (u0, v0) = (l[0].add_r(&-&c[0]), l[1].add_r(&-&c[1]));
            let a = &ld[0] * &ld[0] + &ld[1] * &ld[1];
            let b = u0.scale(&ld[0]).add(&v0.scale(&ld[1])).scale(&int(2));
            let cc = u0.mul(&u0).add(&v0.mul(&v0)).add_r(&-(r * r));
            if a == zero() {
                return Ok(if cc.sign() == Ordering::Equal {
                    EdgeMeet::Along
                } else {
                    EdgeMeet::None
                });
            }
            if !b.is_rational() || !cc.is_rational() {
                // A generatrix's base is irrational only for sections, not
                // for edges.
                return Err(Error::ComputationLimit(
                    "an irrational line against a cylinder",
                ));
            }
            let (Some(b), Some(cc)) = (b.rational(), cc.rational()) else {
                unreachable!("rational above")
            };
            let roots = quadratic(&a, b, cc)?;
            Ok(EdgeMeet::Points(
                roots
                    .into_iter()
                    .map(|t| {
                        let x = qadd(p, &qscale(d, &t));
                        (Pos::T(t), x)
                    })
                    .collect(),
            ))
        }
        (Crv::Conic { c, a, b }, Surf::Plane { p: p0, m }) => {
            let alpha = dot(m, a);
            let beta = dot(m, b);
            let gamma = dot(m, &sub(p0, c));
            match trig(&alpha, &beta, &gamma)? {
                None => Ok(EdgeMeet::Along),
                Some(roots) => Ok(EdgeMeet::Points(
                    roots
                        .into_iter()
                        .map(|cs| {
                            let x = conic_point(c, a, b, &cs);
                            (Pos::Ang(cs), x)
                        })
                        .collect(),
                )),
            }
        }
        (Crv::Conic { c, a, b }, Surf::Cyl { c: cy, r: ry, .. }) if own.is_none() => {
            // A cone's rim against a cylinder (S9d.3b): algebraic points.
            Ok(EdgeMeet::Points(super::algebraic::circle_points(
                c, a, b, &py.f, cy, ry,
            )?))
        }
        (Crv::Conic { .. }, Surf::Cyl { c: cy, r: ry, .. }) => {
            // An arc of the edge's own cylinder (`own`) against the other's.
            let (px, cx, rx) = own.expect("an arc edge's own cylinder");
            let Crv::Conic { c, a, b } = curve else {
                unreachable!("matched above")
            };
            match pair.expect("a cylinder pair's relation") {
                CylPair::Mixed(_) | CylPair::Plane { .. } => {
                    unreachable!("a cylinder and a sphere's or two cones'")
                }
                CylPair::Apart => Ok(EdgeMeet::None),
                CylPair::Same => Ok(EdgeMeet::Along),
                CylPair::Parallel { c2, r2 } => {
                    let cs = circle_circle(cx, rx, c2, r2)?;
                    let _ = px;
                    Ok(EdgeMeet::Points(
                        cs.into_iter()
                            .map(|e| {
                                let x = conic_point(c, a, b, &e);
                                (Pos::Ang(e), x)
                            })
                            .collect(),
                    ))
                }
                // A quartic's or parallel lines' points on the circle:
                // algebraic (S9c.2b.2).
                CylPair::Quartic(_) | CylPair::Lines(..) => Ok(EdgeMeet::Points(
                    super::algebraic::circle_points(c, a, b, &py.f, cy, ry)?,
                )),
                CylPair::Crossing(x) => {
                    let mut out = Vec::new();
                    for (p0, m) in &x.planes {
                        let alpha = dot(m, a);
                        let beta = dot(m, b);
                        let gamma = dot(m, &sub(p0, c));
                        if let Some(roots) = trig(&alpha, &beta, &gamma)? {
                            for cs in roots {
                                let x = conic_point(c, a, b, &cs);
                                out.push((Pos::Ang(cs), x));
                            }
                        }
                    }
                    Ok(EdgeMeet::Points(out))
                }
            }
        }
    }
}

/// The unit-free tangent of a curve at a place (the point `x`).
pub(super) fn tangent(curve: &Crv, pos: &Pos, x: &QV) -> QV {
    match (curve, pos) {
        (Crv::Line { d, .. }, _) => qv(d),
        (Crv::Circle(c), _) => c.tangent(x),
        (Crv::Rise(c), _) => c.tangent(x),
        (Crv::Conic { a, b, .. }, Pos::Ang(cs)) => conic_tangent(a, b, cs),
        (Crv::Meet(m), Pos::Ang(cs)) => m.tangent(cs, x),
        (Crv::Cone(c), Pos::Ang(cs)) => c.tangent(cs),
        (Crv::Torus(c), _) => c.tangent(x),
        (Crv::Conic { .. } | Crv::Meet(_) | Crv::Cone(_), Pos::T(_)) => {
            unreachable!("a conic's place is an angle")
        }
    }
}

/// A conic's `(cos, sin)` of a point on it: from `c + a cos + b sin = x`,
/// solved in the plane of `a` and `b`.
pub(super) fn conic_angle(c: &V, a: &V, b: &V, x: &QV) -> [Qd; 2] {
    let d = qsub(x, &qv(c));
    let (aa, ab, bb) = (dot(a, a), dot(a, b), dot(b, b));
    let (da, db) = (qdot(&d, a), qdot(&d, b));
    let det = &aa * &bb - &ab * &ab;
    let inv = int(1) / det;
    [
        da.scale(&bb).sub(&db.scale(&ab)).scale(&inv),
        db.scale(&aa).sub(&da.scale(&ab)).scale(&inv),
    ]
}

/// Whether the circle `c + a cos + b sin` lies beyond a cone's end planes
/// (its heights `w0 +- sqrt(alpha^2 + beta^2)` all above `h` or below 0).
fn beyond_ends(f: &Affine, h: &R, c: &V, a: &V, b: &V) -> bool {
    let row = f.row(2);
    let w0 = dot(row, &sub(c, &f.o));
    let reach = {
        let (alpha, beta) = (dot(row, a), dot(row, b));
        &alpha * &alpha + &beta * &beta
    };
    let above = &w0 - h;
    (above > zero() && &above * &above > reach) || (w0 < zero() && &w0 * &w0 > reach)
}
