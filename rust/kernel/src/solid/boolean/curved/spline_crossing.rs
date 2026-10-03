//! S9f.2b.1: spline walls against cylinder walls on crossing axes
//! (REVIEW_NOTES.md, "S9f.2b refined, before its code").
//!
//! Along a spline wall's ruling at the run parameter `tau`, `X = o + S_x(tau)
//! x + S_y(tau) y + w n` (the spline prism's exact model), the cylinder's
//! function in its own frame's exact rows (`procedural::Other`: `sum_i (g_i
//! . X - e_i)^2 - r^2`, an elliptic cylinder in the world where the stored
//! axes are not exactly orthonormal) is `F = A w^2 + 2 B(tau) w + C(tau)`:
//! with `P_i = g_i . (o + x S_x + y S_y) - e_i` and `q_i = g_i . n`, `A = sum
//! q_i^2` (positive: the axes cross), `B = sum q_i P_i` (degree `p` on each
//! Bézier arc) and `C = sum P_i^2 - r^2` (degree `2 p`). The meeting's
//! branches are `w = (-B +- sqrt(D)) / A`, `D = B^2 - A C` of degree `2 p`,
//! its roots the turning points (the ruling tangent to the cylinder).
//!
//! The section of a spline wall and a cylinder is each branch over each
//! maximal range of the segment's run where `D > 0` (`Crv::WallMeet`, a graph
//! over `tau`): a turning point ends a range where it lies outside either
//! face (the pieces beside it are outside too); one inside both faces is
//! S9f.2b.2's (loops, `OutOfDomain`); one on a face's boundary, at an
//! interior knot or a segment's end, or of multiplicity above one (the
//! cylinder tangent to the wall) is `Degenerate`. A point is on a branch
//! when its profile point lies on the segment within the range, it lies on
//! the cylinder and `A w + B(tau) = sum q_i (g_i . X - e_i)` (half `F`'s
//! derivative in `w`) has the branch's sign, exactly; a piece's midpoint at
//! a rational `tau` lies in `Q(sqrt(D(tau)))`.
//!
//! Vertices: a curve over a spline segment (a cap edge or a crease, `w = h0
//! + h1 S_x + h2 S_y`) meets the cylinder at the roots of `F` along it,
//! degree `2 p` (`wallcrv_cyl`); a cylinder's cap circle meets the wall at
//! the roots of `F` along its cap plane's crease on the wall, degree `2 p`,
//! at its angle there (`conic_wall`; a cap plane holding the wall's axis
//! direction meets the wall in generatrices whose points with the circle lie
//! in `Q(alpha)(sqrt(delta))`, a tower: S9f.2b.2's); the cylinder prism's
//! vertical edges meet the wall by S9f.1's `line_wall` (degree `p`), the
//! spline prism's the cylinder at quadratic surds. Every root is the arc's
//! own parameter (S9f.1's generators), so every vertex on the wall lies in
//! `Q(alpha)` of degree at most `2 p` or in `Q(sqrt(d))`.
use super::meet::{conic_angle, EdgeMeet, Pos, Section};
use super::model::{Affine, Crv, FaceKind, Loc, Prism, Seg, P2};
use super::num::*;
use super::procedural::{other_of, Other};
use super::spline_walls::{combine, qpoint, trim, BArc, Roots, SplineSeg, WallCrv};
use crate::solid::split::{rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

fn loops() -> Error {
    Error::OutOfDomain(
        "a spline wall's meeting with a cylinder turning back inside the faces (S9f.2b.2)",
    )
}

fn tower() -> Error {
    Error::OutOfDomain(
        "a cylinder's cap circle meeting a spline wall in a plane along the wall's axis (S9f.2b.2)",
    )
}

fn tangent_wall() -> Error {
    Error::Degenerate("a cylinder tangent to a spline wall")
}

fn turning_knot() -> Error {
    Error::Degenerate("a spline wall's meeting with a cylinder turning back at a knot")
}

fn turning_edge() -> Error {
    Error::Degenerate("a spline wall's meeting with a cylinder turning back on a face's boundary")
}

fn tangent_curve() -> Error {
    Error::Degenerate("a spline edge of one input tangent to a face of the other")
}

// ------------------------------------------------------------ polynomials

fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if *x == zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

fn padd(a: &[R], b: &[R]) -> Vec<R> {
    let mut out = vec![zero(); a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        out[i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] += x;
    }
    trim(out)
}

fn pscale(a: &[R], k: &R) -> Vec<R> {
    trim(a.iter().map(|x| x * k).collect())
}

/// The cylinder's rows along a curve over an arc, `w = h0 + h1 S_x + h2 S_y`
/// (`h = None`: free, the ruling's `w` coefficient kept apart): per row `P_i
/// + q_i w`, `(P_i(s), q_i)`.
fn rows(f: &Affine, arc: &BArc, other: &Other, h: Option<&[R; 3]>) -> Vec<(Vec<R>, R)> {
    other
        .g
        .iter()
        .zip(&other.e)
        .map(|(g, e)| {
            let q = dot(g, &f.n);
            let (gx, gy) = (dot(g, &f.x), dot(g, &f.y));
            let c0 = dot(g, &f.o) - e;
            match h {
                None => (combine(&c0, &gx, &arc.x, &gy, &arc.y), q),
                Some(h) => (
                    combine(
                        &(&c0 + &q * &h[0]),
                        &(&gx + &q * &h[1]),
                        &arc.x,
                        &(&gy + &q * &h[2]),
                        &arc.y,
                    ),
                    zero(),
                ),
            }
        })
        .collect()
}

/// `A`, `B(s)` and `C(s)` of the cylinder's function along the arc's
/// rulings.
fn ruling_coeffs(f: &Affine, arc: &BArc, other: &Other) -> (R, Vec<R>, Vec<R>) {
    let rs = rows(f, arc, other, None);
    let a = rs.iter().fold(zero(), |acc, (_, q)| acc + q * q);
    let b = rs
        .iter()
        .fold(Vec::new(), |acc: Vec<R>, (p, q)| padd(&acc, &pscale(p, q)));
    let c = rs
        .iter()
        .fold(vec![-(&other.r * &other.r)], |acc: Vec<R>, (p, _)| {
            padd(&acc, &pmul(p, p))
        });
    (a, b, c)
}

/// The cylinder's function along a curve over an arc (`w` given by `h`).
fn along_curve(f: &Affine, arc: &BArc, other: &Other, h: &[R; 3]) -> Vec<R> {
    rows(f, arc, other, Some(h))
        .iter()
        .fold(vec![-(&other.r * &other.r)], |acc: Vec<R>, (p, _)| {
            padd(&acc, &pmul(p, p))
        })
}

/// `A`, `B` and `C` at a rational profile point.
fn coeffs_at(f: &Affine, other: &Other, uv: &[R; 2]) -> (R, R, R) {
    let base = add(&add(&f.o, &scale(&f.x, &uv[0])), &scale(&f.y, &uv[1]));
    let (mut a, mut b, mut c) = (zero(), zero(), -(&other.r * &other.r));
    for (g, e) in other.g.iter().zip(&other.e) {
        let (p, q) = (dot(g, &base) - e, dot(g, &f.n));
        a += &q * &q;
        b += &q * &p;
        c += &p * &p;
    }
    (a, b, c)
}

// ------------------------------------------------------------ the curve

/// A branch of a spline wall's meeting with a cylinder over a range of the
/// segment's run (S9f.2b.1), placed by the run parameter.
#[derive(Debug, Clone)]
pub(super) struct WallMeetCrv {
    pub(super) seg: Arc<SplineSeg>,
    /// The spline prism's frame (its model's).
    pub(super) f: Affine,
    pub(super) other: Other,
    /// The spline prism's operand: its face in the section is the wall.
    pub(super) carrier: usize,
    /// The branch: `A w + B(tau)` positive or negative.
    pub(super) plus: bool,
    pub(super) range: [Qd; 2],
}

impl PartialEq for WallMeetCrv {
    fn eq(&self, o: &Self) -> bool {
        // One wall, one cylinder (a cylinder's faces split at seams share
        // it), one branch, one range.
        Arc::ptr_eq(&self.seg, &o.seg)
            && self.carrier == o.carrier
            && self.other.g == o.other.g
            && self.other.e == o.other.e
            && self.other.r == o.other.r
            && self.plus == o.plus
            && self.range[0].cmp(&o.range[0]) == Ordering::Equal
            && self.range[1].cmp(&o.range[1]) == Ordering::Equal
    }
}

impl WallMeetCrv {
    /// Half `F`'s derivative in `w` at a point: `sum q_i (g_i . X - e_i)`.
    fn slope(&self, x: &QV) -> Qd {
        self.other
            .g
            .iter()
            .zip(&self.other.e)
            .fold(Qd::rat(zero()), |acc, (g, e)| {
                acc.add(&qdot(x, g).add_r(&-e.clone()).scale(&dot(g, &self.f.n)))
            })
    }

    fn within(&self, tau: &Qd) -> bool {
        tau.cmp(&self.range[0]) != Ordering::Less && tau.cmp(&self.range[1]) != Ordering::Greater
    }

    /// The run parameter of a point on the curve, if it is on it.
    pub(super) fn locate(&self, x: &QV) -> Option<Qd> {
        let l = self.f.local_q(x);
        let tau = self.seg.locate(&[l[0].clone(), l[1].clone()])?;
        if !self.within(&tau) || self.other.value(x).sign() != Ordering::Equal {
            return None;
        }
        let want = if self.plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        (self.slope(x).sign() == want).then_some(tau)
    }

    pub(super) fn on(&self, x: &QV) -> bool {
        self.locate(x).is_some()
    }

    pub(super) fn place(&self, x: &QV) -> Qd {
        self.locate(x)
            .expect("a point on a spline wall's meeting at a known parameter")
    }

    /// The curve's point at a rational run parameter inside its range.
    pub(super) fn at(&self, tau: &R) -> Option<QV> {
        let uv = self.seg.point(&Qd::rat(tau.clone()));
        let (Some(u), Some(v)) = (uv[0].rational(), uv[1].rational()) else {
            return None;
        };
        let uv = [u.clone(), v.clone()];
        let (a, b, c) = coeffs_at(&self.f, &self.other, &uv);
        let d = &b * &b - &a * &c;
        if d <= zero() {
            return None;
        }
        let inv = int(1) / &a;
        let s = if self.plus { int(1) } else { int(-1) };
        let w = Qd::new(-&b * &inv, s * &inv, d);
        Some(qpoint(
            &self.f,
            &Qd::rat(uv[0].clone()),
            &Qd::rat(uv[1].clone()),
            &w,
        ))
    }

    /// The tangent along increasing `tau` at a point (unit-free): `F_w dS -
    /// F_tau n` on the plus branch (`F_w > 0`), its opposite on the other.
    pub(super) fn tangent(&self, x: &QV, tau: &Qd) -> QV {
        let d = self.seg.deriv(tau);
        let ds = qadd(&qscale(&self.f.x, &d[0]), &qscale(&self.f.y, &d[1]));
        let mut ft = Qd::rat(zero());
        for (g, e) in self.other.g.iter().zip(&self.other.e) {
            ft = ft.add(&qdot(x, g).add_r(&-e.clone()).mul(&qdot(&ds, g)));
        }
        let fw = self.slope(x);
        let t: QV = std::array::from_fn(|k| ds[k].mul(&fw).sub(&ft.scale(&self.f.n[k])));
        if self.plus {
            t
        } else {
            t.map(|c| c.neg())
        }
    }

    /// Points in binary64 from `t0` to `t1` (run parameters).
    pub(super) fn samples(&self, t0: f64, t1: f64, n: usize) -> Vec<[f64; 3]> {
        let fl = |v: &V| v.clone().map(|x| rational_f64(&x));
        let (o, x, y, nn) = (fl(&self.f.o), fl(&self.f.x), fl(&self.f.y), fl(&self.f.n));
        let rows: Vec<([f64; 3], f64)> = self
            .other
            .g
            .iter()
            .zip(&self.other.e)
            .map(|(g, e)| (fl(g), rational_f64(e)))
            .collect();
        let r = rational_f64(&self.other.r);
        (0..=n)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let [u, v] = self.seg.point_f64(t);
                let base: [f64; 3] = [0, 1, 2].map(|k| o[k] + x[k] * u + y[k] * v);
                let (mut a, mut b, mut c) = (0.0, 0.0, -r * r);
                for (g, e) in &rows {
                    let p = g[0] * base[0] + g[1] * base[1] + g[2] * base[2] - e;
                    let q = g[0] * nn[0] + g[1] * nn[1] + g[2] * nn[2];
                    a += q * q;
                    b += q * p;
                    c += p * p;
                }
                let s = if self.plus { 1.0 } else { -1.0 };
                let sq = s * (b * b - a * c).max(0.0).sqrt();
                let (pp, mm) = (-b + sq, -b - sq);
                let w = if pp.abs() >= mm.abs() { pp / a } else { c / mm };
                [0, 1, 2].map(|k| base[k] + nn[k] * w)
            })
            .collect()
    }
}

// ------------------------------------------------------------ the section

/// A spline wall's section by a cylinder wall of another prism on a
/// crossing axis (`sm`'s face `sf` holding segment `seg`, `cm`'s face `cf`
/// on the cylinder `c`, `r`): each branch over each range of the run where
/// the discriminant is positive.
#[allow(clippy::too_many_arguments)]
pub(super) fn section(
    sm: &Prism,
    sf: usize,
    seg: &Arc<SplineSeg>,
    cm: &Prism,
    cf: usize,
    c: &P2,
    r: &R,
) -> Result<Section> {
    let other = other_of(&cm.f, c, r);
    let f = &sm.f;
    let rs = match seg.roots_of(&|arc: &BArc| {
        let (a, b, cc) = ruling_coeffs(f, arc, &other);
        padd(&pmul(&b, &b), &pscale(&cc, &-a))
    })? {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(tangent_wall()),
    };
    let first = Qd::rat(seg.first.clone());
    let last = Qd::rat(seg.last.clone());
    let mut bounds = vec![first.clone()];
    for root in &rs {
        // The turning point: on the ruling at w = -B / A.
        let uv = seg.point(&root.tau);
        let base = qpoint(f, &uv[0], &uv[1], &Qd::rat(zero()));
        let (mut a, mut b) = (zero(), Qd::rat(zero()));
        for (g, e) in other.g.iter().zip(&other.e) {
            let q = dot(g, &f.n);
            a += &q * &q;
            b = b.add(&qdot(&base, g).add_r(&-e.clone()).scale(&q));
        }
        let w = b.scale(&(int(-1) / a));
        let x = qpoint(f, &uv[0], &uv[1], &w);
        let inside = (sm.in_face(sf, &x), cm.in_face(cf, &x));
        let near = !matches!(inside, (Loc::Out, _) | (_, Loc::Out));
        // Outside a face but within the resolution of both: the meeting's
        // end there a rounding away from turning back, its square root's
        // slope unbounded (a frame normalized to an ulp off a tangent
        // ruling on a cap's edge).
        let tol = sm.tolerance.linear().max(cm.tolerance.linear());
        let close = |pr: &Prism, fi: usize, at: Loc| at != Loc::Out || outside(pr, fi, &x) <= tol;
        if !near && close(sm, sf, inside.0) && close(cm, cf, inside.1) {
            return Err(turning_edge());
        }
        if near {
            if root.mult > 1 {
                return Err(tangent_wall());
            }
            if root.knot {
                return Err(turning_knot());
            }
            if root.end || inside != (Loc::In, Loc::In) {
                return Err(turning_edge());
            }
            return Err(loops());
        }
        if !root.end {
            bounds.push(root.tau.clone());
        }
    }
    bounds.push(last);
    let mut out = Vec::new();
    for w in bounds.windows(2) {
        let k = super::graph::rational_between_num(&w[0], &w[1])?;
        let uv = seg.point(&Qd::rat(k));
        let (Some(u), Some(v)) = (uv[0].rational(), uv[1].rational()) else {
            unreachable!("a rational run parameter's point is rational")
        };
        let (a, b, cc) = coeffs_at(f, &other, &[u.clone(), v.clone()]);
        if &b * &b - &a * &cc <= zero() {
            continue;
        }
        for plus in [true, false] {
            out.push(Crv::WallMeet(Box::new(WallMeetCrv {
                seg: seg.clone(),
                f: f.clone(),
                other: other.clone(),
                carrier: seg.op,
                plus,
                range: [w[0].clone(), w[1].clone()],
            })));
        }
    }
    Ok(Section::Curves(out))
}

/// How far (binary64, in the prism's local units) a point on a wall face's
/// surface lies outside the face: past its heights and, where its generatrix
/// leaves the face's line or arc, from the nearer of their ends.
fn outside(pr: &Prism, fi: usize, x: &QV) -> f64 {
    let FaceKind::Wall(b, j) = pr.faces[fi].kind else {
        return f64::INFINITY;
    };
    let l = pr.f.local_q(x);
    let w = if l[2].cmp(&Qd::rat(pr.lo.clone())) == Ordering::Less {
        Qd::rat(pr.lo.clone())
    } else if l[2].cmp(&Qd::rat(pr.hi.clone())) == Ordering::Greater {
        Qd::rat(pr.hi.clone())
    } else {
        l[2].clone()
    };
    let height = (l[2].to_f64() - w.to_f64()).abs();
    let across = if pr.in_face(fi, &qpoint(&pr.f, &l[0], &l[1], &w)) == Loc::Out {
        let at = [l[0].to_f64(), l[1].to_f64()];
        let gap = |e: &P2| (at[0] - rational_f64(&e[0])).hypot(at[1] - rational_f64(&e[1]));
        let seg = &pr.bounds[b].segs[j];
        if matches!(seg, Seg::Spline(_)) {
            0.0
        } else {
            gap(seg.start()).min(gap(seg.end()))
        }
    } else {
        0.0
    };
    height.hypot(across)
}

// ------------------------------------------------------------ the vertices

/// Where a curve over a spline segment (a cap edge or a crease) meets a
/// cylinder wall of another prism on a crossing axis (`c`, `r` on frame
/// `cf`): the roots of the cylinder's function along it.
pub(super) fn wallcrv_cyl(curve: &WallCrv, cf: &Affine, c: &P2, r: &R) -> Result<EdgeMeet> {
    let other = other_of(cf, c, r);
    let rs = match curve
        .seg
        .roots_of(&|arc: &BArc| along_curve(&curve.f, arc, &other, &curve.h))?
    {
        Roots::At(rs) => rs,
        Roots::Along => return Ok(EdgeMeet::Along),
        Roots::Partly => return Err(tangent_curve()),
    };
    let mut out = Vec::new();
    for root in rs {
        if root.mult > 1 {
            return Err(tangent_curve());
        }
        let x = curve.point(&root.tau);
        out.push((Pos::T(root.tau), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// Where a cylinder's cap edge (the conic `cc + a cos + b sin` on its own
/// cylinder `c`, `r` of the model `cm`) meets a spline wall (segment `seg` on
/// frame `sf`) on a crossing axis: the cylinder's function along the cap
/// plane's crease on the wall, each root at its angle on the circle.
#[allow(clippy::too_many_arguments)]
pub(super) fn conic_wall(
    cm: &Prism,
    c: &P2,
    r: &R,
    cc: &V,
    a: &V,
    b: &V,
    sf: &Affine,
    seg: &SplineSeg,
) -> Result<EdgeMeet> {
    let m = cross(a, b);
    let mn = dot(&m, &sf.n);
    let k = dot(&m, &sub(&sf.o, cc));
    let (mx, my) = (dot(&m, &sf.x), dot(&m, &sf.y));
    if mn == zero() {
        // The cap plane holds the wall's axis direction: its generatrices
        // on the wall, the circle's points there in a tower (S9f.2b.2's);
        // none when the plane misses the segment.
        return match seg.roots_of(&|arc: &BArc| combine(&k, &mx, &arc.x, &my, &arc.y))? {
            Roots::At(rs) if rs.is_empty() => Ok(EdgeMeet::None),
            _ => Err(tower()),
        };
    }
    let inv = int(-1) / &mn;
    let h = [&k * &inv, &mx * &inv, &my * &inv];
    let other = other_of(&cm.f, c, r);
    let rs = match seg.roots_of(&|arc: &BArc| along_curve(sf, arc, &other, &h))? {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(tangent_curve()),
    };
    let mut out = Vec::new();
    for root in rs {
        if root.mult > 1 {
            return Err(tangent_curve());
        }
        let uv = seg.point(&root.tau);
        let w = uv[0].scale(&h[1]).add(&uv[1].scale(&h[2])).add_r(&h[0]);
        let x = qpoint(sf, &uv[0], &uv[1], &w);
        let cs = conic_angle(cc, a, b, &x);
        out.push((Pos::Ang(cs), x));
    }
    Ok(EdgeMeet::Points(out))
}
