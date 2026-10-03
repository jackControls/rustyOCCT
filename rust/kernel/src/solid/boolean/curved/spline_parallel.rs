//! S9f.2a: spline walls against arc, circle and spline walls of a prism
//! whose axis is exactly parallel (REVIEW_NOTES.md, "S9f.2a refined, before
//! its code").
//!
//! Both prisms' walls run along one direction, so every meeting of a spline
//! wall with the other's curved wall is generatrices, over the 2D crossings
//! of the profiles' curves projected along the axis. The other frame's
//! `(u', v')` are rational affine functions of this frame's `(u, v)` (the
//! axis carries none of them: `map2`), so
//!
//! * a spline arc `S(s)` (degree `p`) meets a cylinder `(u' - c_u)^2 + (v' -
//!   c_v)^2 = r^2` (a circle in its own frame, an ellipse in this one when
//!   the map is no similarity) at the roots of `E(s) = |m(S(s)) - c|^2 -
//!   r^2`, degree `2 p`, exact;
//! * two spline arcs meet where the first's parameter is a root of `R(s) =
//!   f(m(S(s)))`, `f` the second arc's implicit equation in its own frame
//!   (`Implicit`: the Sylvester resultant `Res_sigma(B_x(sigma) - X,
//!   B_y(sigma) - Y)`, S9a.2's), degree at most `p q`; the second's
//!   parameter there is `sigma = -b / a` of the first subresultant `a sigma
//!   + b` of `B_x(sigma) - X` and `B_y(sigma) - Y` at the point (`a`, `b`
//!   polynomials in `(X, Y)`, both from exact evaluations on a grid and
//!   interpolation), an element of the root's own field `Q(alpha)`; a
//!   crossing is kept when `sigma` lies in `[0, 1]`, and `a = 0` there (a
//!   node or cusp of the second arc's curve) is `ComputationLimit`.
//!
//! Every root is the spline arc's own parameter (S9f.1's generators): a
//! crossing with a cylinder in the spline's field, two splines' in the
//! object's (operand A's) spline's, whichever way the crossing is found, so
//! a point found twice is one number of one field. The cap edges and
//! creases over a spline segment meet the other's curved wall at the same
//! roots (they project onto the segment); a cylinder's cap edge (a conic)
//! meets a spline wall at the same roots too, at its angle on its circle.
//! Tangencies are `Degenerate`: a root of multiplicity above one, or a
//! crossing at a knot of either curve where the legs beside the knot do
//! not lie on the other curve's two sides; arcs whose degrees' product
//! exceeds 16 are `ComputationLimit`, and a spline arc on the other's curve
//! (a resultant identically zero) is refused as coincident walls.
use super::meet::{conic_point, EdgeMeet, Pos, Section};
use super::model::{Affine, Crv, Prism, P2};
use super::num::*;
use super::spline_walls::{combine, peval, peval_r, qpoint, trim, Root, Roots, SplineSeg, WallCrv};
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

/// The largest product of two spline arcs' degrees whose crossings are
/// computed: their field's degree (REVIEW_NOTES.md, S9f.2a's (4)).
const DEGREE_PRODUCT: usize = 16;

fn crossing_axes_cylinder() -> Error {
    Error::OutOfDomain("a spline wall against a cylinder on crossing axes (S9f.2b)")
}

fn crossing_axes_splines() -> Error {
    Error::OutOfDomain("spline walls against spline walls on crossing axes (refused, S9f)")
}

fn coincident() -> Error {
    Error::OutOfDomain(
        "spline walls of both inputs on one surface in different frames (refused, S9f)",
    )
}

fn tangent_cylinder() -> Error {
    Error::Degenerate("a spline wall tangent to a cylinder along a generatrix")
}

fn tangent_splines() -> Error {
    Error::Degenerate("two spline walls tangent along a generatrix")
}

fn touching_knot() -> Error {
    Error::Degenerate("a curved wall touching a spline wall at a knot's generatrix")
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

fn degree(p: &[R]) -> usize {
    trim(p.to_vec()).len().saturating_sub(1)
}

/// The determinant of a square rational matrix (Gaussian elimination).
fn det(mut m: Vec<Vec<R>>) -> R {
    let n = m.len();
    let mut out = int(1);
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| m[r][col] != zero()) else {
            return zero();
        };
        if pivot != col {
            m.swap(pivot, col);
            out = -out;
        }
        let row = m[col].clone();
        for r in m.iter_mut().skip(col + 1) {
            if r[col] == zero() {
                continue;
            }
            let f = &r[col] / &row[col];
            for (x, p) in r.iter_mut().zip(&row).skip(col) {
                *x -= &f * p;
            }
        }
        out *= &row[col];
    }
    out
}

/// The coefficients (ascending) of the polynomial through `(i, ys[i])`, `i
/// = 0, 1, ...` (Newton's divided differences).
fn interpolate(ys: &[R]) -> Vec<R> {
    let n = ys.len();
    let xs: Vec<R> = (0..n).map(|i| int(i as i64)).collect();
    let mut coef = ys.to_vec();
    for j in 1..n {
        for i in (j..n).rev() {
            coef[i] = (&coef[i] - &coef[i - 1]) / (&xs[i] - &xs[i - j]);
        }
    }
    let mut out = vec![zero(); n];
    for i in (0..n).rev() {
        let mut next = vec![zero(); n];
        for (k, c) in out.iter().enumerate() {
            if *c == zero() {
                continue;
            }
            if k + 1 < n {
                next[k + 1] += c;
            }
            next[k] -= c * &xs[i];
        }
        next[0] += &coef[i];
        out = next;
    }
    out
}

/// A polynomial in two variables, `c[i][j] X^i Y^j`.
#[derive(Debug, Clone)]
pub(super) struct Bivar(Vec<Vec<R>>);

impl Bivar {
    /// The polynomial taking `value(i, j)` at the integers `0 <= i < nx`,
    /// `0 <= j < ny` (of degree below `nx` in `X`, `ny` in `Y`).
    fn interpolate(nx: usize, ny: usize, value: impl Fn(&R, &R) -> R) -> Self {
        // In Y for each X, then in X for each power of Y.
        let rows: Vec<Vec<R>> = (0..nx)
            .map(|i| {
                let x = int(i as i64);
                let ys: Vec<R> = (0..ny).map(|j| value(&x, &int(j as i64))).collect();
                let mut c = interpolate(&ys);
                c.resize(ny, zero());
                c
            })
            .collect();
        let mut out = vec![vec![zero(); ny]; nx];
        for j in 0..ny {
            let ys: Vec<R> = rows.iter().map(|r| r[j].clone()).collect();
            for (i, c) in interpolate(&ys).into_iter().enumerate() {
                out[i][j] = c;
            }
        }
        Self(out)
    }

    /// Its value at a point of any one field.
    pub(super) fn eval(&self, x: &Qd, y: &Qd) -> Qd {
        self.0.iter().rev().fold(Qd::rat(zero()), |acc, row| {
            acc.mul(x).add(
                &row.iter()
                    .rev()
                    .fold(Qd::rat(zero()), |a, c| a.mul(y).add_r(c)),
            )
        })
    }

    /// `f(X(t), Y(t))` for polynomials `X`, `Y` in `t`.
    fn compose(&self, x: &[R], y: &[R]) -> Vec<R> {
        let mut out: Vec<R> = Vec::new();
        for row in self.0.iter().rev() {
            let inner = row.iter().rev().fold(Vec::new(), |a: Vec<R>, c| {
                padd(&pmul(&a, y), std::slice::from_ref(c))
            });
            out = padd(&pmul(&out, x), &inner);
        }
        out
    }
}

/// The Sylvester resultant of two polynomials in one variable (ascending,
/// their degrees `m, n >= 1`).
fn resultant(f: &[R], g: &[R]) -> R {
    let (m, n) = (f.len() - 1, g.len() - 1);
    let size = m + n;
    let mut rows = vec![vec![zero(); size]; size];
    for (r, row) in rows.iter_mut().enumerate().take(n) {
        for k in 0..=m {
            row[r + k] = f[m - k].clone();
        }
    }
    for (r, row) in rows.iter_mut().skip(n).enumerate() {
        for k in 0..=n {
            row[r + k] = g[n - k].clone();
        }
    }
    det(rows)
}

/// The first subresultant `a x + b` of two polynomials of degrees `m, n >=
/// 1`, `m + n >= 3` (ascending): `(a, b)`.
fn subresultant1(f: &[R], g: &[R]) -> (R, R) {
    let (m, n) = (f.len() - 1, g.len() - 1);
    let rows_n = m + n - 2;
    let cols = m + n - 1;
    // Columns by the powers m + n - 2 down to 0.
    let mut rows = Vec::with_capacity(rows_n);
    for s in (0..n - 1).rev() {
        let mut row = vec![zero(); cols];
        for (e, c) in f.iter().enumerate() {
            row[cols - 1 - (e + s)] = c.clone();
        }
        rows.push(row);
    }
    for s in (0..m - 1).rev() {
        let mut row = vec![zero(); cols];
        for (e, c) in g.iter().enumerate() {
            row[cols - 1 - (e + s)] = c.clone();
        }
        rows.push(row);
    }
    let pick = |last: usize| -> Vec<Vec<R>> {
        rows.iter()
            .map(|r| {
                let mut out: Vec<R> = r[..cols - 2].to_vec();
                out.push(r[last].clone());
                out
            })
            .collect()
    };
    (det(pick(cols - 2)), det(pick(cols - 1)))
}

/// How an arc's parameter is read from a point on its curve.
#[derive(Debug, Clone)]
enum Inversion {
    /// `sigma = -b / a` (the first subresultant).
    Sub(Bivar, Bivar),
    /// A coordinate of degree one: `sigma = (coordinate - c0) / c1`.
    Linear(usize, R, R),
    /// None (a straight arc run along an axis by a nonlinear parameter).
    Unknown,
}

/// An arc's implicit equation in its own frame and its parameter's
/// inversion (S9f.2a).
#[derive(Debug, Clone)]
pub(super) struct Implicit {
    f: Bivar,
    inversion: Inversion,
}

impl Implicit {
    /// Of the arc `(x(s), y(s))` (power coefficients, ascending).
    pub(super) fn of(x: &[R], y: &[R]) -> Self {
        let (x, y) = (trim(x.to_vec()), trim(y.to_vec()));
        let pad = |p: Vec<R>| if p.is_empty() { vec![zero()] } else { p };
        let (x, y) = (pad(x), pad(y));
        let (dx, dy) = (degree(&x), degree(&y));
        let shifted = |p: &[R], v: &R| {
            let mut q = p.to_vec();
            q[0] -= v;
            q
        };
        let f = if dx == 0 {
            // A vertical straight arc: its line `X = x0`.
            Bivar(vec![vec![-&x[0]], vec![int(1)]])
        } else if dy == 0 {
            Bivar(vec![vec![-&y[0], int(1)]])
        } else {
            Bivar::interpolate(dy + 1, dx + 1, |u, v| {
                resultant(&shifted(&x, u), &shifted(&y, v))
            })
        };
        let inversion = if dx == 1 {
            Inversion::Linear(0, x[0].clone(), x[1].clone())
        } else if dy == 1 {
            Inversion::Linear(1, y[0].clone(), y[1].clone())
        } else if dx == 0 || dy == 0 {
            Inversion::Unknown
        } else {
            let a = Bivar::interpolate(dy, dx, |u, v| {
                subresultant1(&shifted(&x, u), &shifted(&y, v)).0
            });
            let b = Bivar::interpolate(dy, dx, |u, v| {
                subresultant1(&shifted(&x, u), &shifted(&y, v)).1
            });
            Inversion::Sub(a, b)
        };
        Self { f, inversion }
    }

    /// The implicit equation's value at a point (of any one field).
    pub(super) fn value(&self, x: &[Qd; 2]) -> Qd {
        self.f.eval(&x[0], &x[1])
    }

    /// The arc's parameter at a point of its curve, or `None` where the
    /// inversion is undetermined (a node or cusp of the curve).
    pub(super) fn param(&self, x: &[Qd; 2]) -> Option<Qd> {
        match &self.inversion {
            Inversion::Linear(i, c0, c1) => Some(x[*i].add_r(&-c0).scale(&(int(1) / c1))),
            Inversion::Sub(a, b) => {
                let av = a.eval(&x[0], &x[1]);
                if av.sign() == Ordering::Equal {
                    return None;
                }
                Some(b.eval(&x[0], &x[1]).neg().mul(&av.recip()?))
            }
            Inversion::Unknown => None,
        }
    }
}

// ------------------------------------------------------------ frames

/// The other frame's `(u', v')` of this frame's `(u, v)`: `[[c, cu, cv];
/// 2]`, `u' = c + cu u + cv v`; `None` unless the axes are exactly
/// parallel.
pub(super) type Map2 = [[R; 3]; 2];

pub(super) fn map2(from: &Affine, to: &Affine) -> Option<Map2> {
    if !is_zero(&cross(&from.n, &to.n)) {
        return None;
    }
    let d = sub(&from.o, &to.o);
    let row = |k: usize| {
        let r = to.row(k);
        [dot(r, &d), dot(r, &from.x), dot(r, &from.y)]
    };
    Some([row(0), row(1)])
}

fn apply(m: &Map2, x: &[Qd; 2]) -> [Qd; 2] {
    let one = |k: usize| {
        x[0].scale(&m[k][1])
            .add(&x[1].scale(&m[k][2]))
            .add_r(&m[k][0])
    };
    [one(0), one(1)]
}

fn apply_dir(m: &Map2, d: &[Qd; 2]) -> [Qd; 2] {
    let one = |k: usize| d[0].scale(&m[k][1]).add(&d[1].scale(&m[k][2]));
    [one(0), one(1)]
}

/// An arc's coordinates mapped: `(u'(s), v'(s))` as polynomials.
fn mapped(m: &Map2, x: &[R], y: &[R]) -> [Vec<R>; 2] {
    [
        combine(&m[0][0], &m[0][1], x, &m[0][2], y),
        combine(&m[1][0], &m[1][1], x, &m[1][2], y),
    ]
}

fn cross2(u: &[Qd; 2], v: &[Qd; 2]) -> Ordering {
    mixed_dot_sign(&[u[0].clone(), u[1].neg()], &[v[1].clone(), v[0].clone()])
}

/// Whether a curve through a knot of a spline crosses it there: the legs
/// arriving at and leaving the knot lie strictly on the two sides of the
/// curve's tangent `t` (all in one frame).
fn crosses_at_knot(arrive: &[Qd; 2], leave: &[Qd; 2], t: &[Qd; 2]) -> bool {
    let before = cross2(t, &[arrive[0].neg(), arrive[1].neg()]);
    let after = cross2(t, leave);
    before != Ordering::Equal && after != Ordering::Equal && before != after
}

/// The legs of a segment at a knot (a run parameter at an arc's end): the
/// derivative at the end of the arc before it and at the start of the arc
/// after it.
fn knot_legs(seg: &SplineSeg, tau: &Qd) -> ([Qd; 2], [Qd; 2]) {
    let k = seg
        .arcs
        .iter()
        .position(|a| tau.cmp(&Qd::rat(a.d[1].clone())) == Ordering::Equal)
        .expect("a knot at an arc's end");
    let (a, b) = (&seg.arcs[k], &seg.arcs[k + 1]);
    let at = |p: &[R], s: i64| Qd::rat(peval_r(p, &int(s)));
    ([at(&a.dx, 1), at(&a.dy, 1)], [at(&b.dx, 0), at(&b.dy, 0)])
}

// ------------------------------------------------------------ a spline against a cylinder

/// The roots along a spline segment (on frame `f`) of a cylinder's equation
/// (`c`, `r` on frame `cf`): `E = |m(S) - c|^2 - r^2` on each arc.
fn cylinder_roots(seg: &SplineSeg, f: &Affine, cf: &Affine, c: &P2, r: &R) -> Result<Vec<Root>> {
    let m = map2(f, cf).ok_or_else(crossing_axes_cylinder)?;
    let rs = seg.roots_of(&|a| {
        let [u, v] = mapped(&m, &a.x, &a.y);
        let (du, dv) = (padd(&u, &[-&c[0]]), padd(&v, &[-&c[1]]));
        padd(&padd(&pmul(&du, &du), &pmul(&dv, &dv)), &[-(r * r)])
    })?;
    let rs = match rs {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(coincident()),
    };
    for root in &rs {
        if root.mult > 1 {
            return Err(tangent_cylinder());
        }
        if root.knot {
            // The circle's tangent at the point, in the spline's frame.
            let uv = seg.point(&root.tau);
            let p = apply(&m, &uv);
            let rel = [p[0].add_r(&-&c[0]), p[1].add_r(&-&c[1])];
            let back = map2(cf, f).expect("parallel both ways");
            let t = apply_dir(&back, &[rel[1].neg(), rel[0].clone()]);
            let (arrive, leave) = knot_legs(seg, &root.tau);
            if !crosses_at_knot(&arrive, &leave, &t) {
                return Err(touching_knot());
            }
        }
    }
    Ok(rs)
}

/// A spline wall's section by a cylinder wall of another prism (`sm`'s
/// segment, `cm`'s cylinder `c`, `r`): generatrices over the crossings
/// strictly inside the segment's run.
pub(super) fn spline_cyl(
    sm: &Prism,
    seg: &Arc<SplineSeg>,
    cm: &Prism,
    c: &P2,
    r: &R,
) -> Result<Section> {
    let rs = cylinder_roots(seg, &sm.f, &cm.f, c, r)?;
    Ok(Section::Curves(
        rs.into_iter()
            .filter(|root| !root.end)
            .map(|root| {
                let uv = seg.point(&root.tau);
                Crv::Line {
                    p: qpoint(&sm.f, &uv[0], &uv[1], &Qd::rat(zero())),
                    d: sm.f.n.clone(),
                }
            })
            .collect(),
    ))
}

/// Where a curve over a spline segment (a cap edge or a crease) meets a
/// cylinder wall of another prism (`c`, `r` on frame `cf`).
pub(super) fn wallcrv_cyl(curve: &WallCrv, cf: &Affine, c: &P2, r: &R) -> Result<EdgeMeet> {
    // S9f.2b: on a crossing axis.
    if map2(&curve.f, cf).is_none() {
        return super::spline_crossing::wallcrv_cyl(curve, cf, c, r);
    }
    let rs = cylinder_roots(&curve.seg, &curve.f, cf, c, r)?;
    Ok(EdgeMeet::Points(
        rs.into_iter()
            .map(|root| {
                let x = curve.point(&root.tau);
                (Pos::T(root.tau), x)
            })
            .collect(),
    ))
}

/// Where a cylinder's cap edge (the conic `cc + a cos + b sin`, on the
/// cylinder `c`, `r` of the model `cm`) meets a spline wall (segment `seg`
/// on frame `sf`): at the cylinder's roots along the segment, each at its
/// angle on the circle.
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
    // S9f.2b: on a crossing axis.
    if map2(sf, &cm.f).is_none() {
        return super::spline_crossing::conic_wall(cm, c, r, cc, a, b, sf, seg);
    }
    let rs = cylinder_roots(seg, sf, &cm.f, c, r)?;
    let m = map2(sf, &cm.f).expect("parallel axes");
    let inv = int(1) / r;
    Ok(EdgeMeet::Points(
        rs.into_iter()
            .map(|root| {
                let p = apply(&m, &seg.point(&root.tau));
                let cs = [
                    p[0].add_r(&-&c[0]).scale(&inv),
                    p[1].add_r(&-&c[1]).scale(&inv),
                ];
                let x = conic_point(cc, a, b, &cs);
                (Pos::Ang(cs), x)
            })
            .collect(),
    ))
}

// ------------------------------------------------------------ two spline walls

/// A crossing of two spline segments' curves (projected along their
/// common axis): each one's run parameter, the point in the first's frame
/// and whether it is at either one's start or end.
#[derive(Debug, Clone)]
pub(super) struct Crossing {
    pub(super) tau: [Qd; 2],
    pub(super) uv: [Qd; 2],
    pub(super) end: [bool; 2],
}

/// The crossings of two spline segments (on frames `fa`, `fb`), in the
/// object's segment's fields whichever is given first.
pub(super) fn crossings(
    sa: &SplineSeg,
    fa: &Affine,
    sb: &SplineSeg,
    fb: &Affine,
) -> Result<Vec<Crossing>> {
    if sa.op <= sb.op {
        return ordered(sa, fa, sb, fb);
    }
    let back = map2(fb, fa).ok_or_else(crossing_axes_splines)?;
    Ok(ordered(sb, fb, sa, fa)?
        .into_iter()
        .map(|c| Crossing {
            uv: apply(&back, &c.uv),
            tau: [c.tau[1].clone(), c.tau[0].clone()],
            end: [c.end[1], c.end[0]],
        })
        .collect())
}

/// The crossings of segment `p` (whose arcs' parameters are the fields)
/// with segment `q`.
fn ordered(p: &SplineSeg, fp: &Affine, q: &SplineSeg, fq: &Affine) -> Result<Vec<Crossing>> {
    let m = map2(fp, fq).ok_or_else(crossing_axes_splines)?;
    let back = map2(fq, fp).expect("parallel both ways");
    let deg = |s: &SplineSeg| s.arcs.iter().map(|a| a.cps.len() - 1).max().unwrap_or(0);
    if deg(p) * deg(q) > DEGREE_PRODUCT {
        return Err(Error::ComputationLimit(
            "two spline walls of degrees whose product exceeds 16 (S9f.2a)",
        ));
    }
    // p's control hull in q's frame, for screening q's arcs.
    let hull: Vec<[R; 2]> = p
        .hull()
        .map(|c| {
            [
                &m[0][0] + &m[0][1] * &c[0] + &m[0][2] * &c[1],
                &m[1][0] + &m[1][1] * &c[0] + &m[1][2] * &c[1],
            ]
        })
        .collect();
    let (one, nil) = (Qd::rat(int(1)), Qd::rat(zero()));
    let n = q.arcs.len();
    let mut out = Vec::new();
    for (k, qa) in q.arcs.iter().enumerate() {
        let apart = |i: usize| {
            let lo = qa.cps.iter().map(|c| &c[i]).min().expect("a pole");
            let hi = qa.cps.iter().map(|c| &c[i]).max().expect("a pole");
            hull.iter().all(|c| &c[i] < lo) || hull.iter().all(|c| &c[i] > hi)
        };
        if apart(0) || apart(1) {
            continue;
        }
        let imp = qa.implicit();
        let rs = match p.roots_of(&|a| {
            let [u, v] = mapped(&m, &a.x, &a.y);
            imp.f.compose(&u, &v)
        })? {
            Roots::At(rs) => rs,
            Roots::Along | Roots::Partly => return Err(coincident()),
        };
        for root in rs {
            let uv = p.point(&root.tau);
            let at = apply(&m, &uv);
            let sigma = imp.param(&at).ok_or(Error::ComputationLimit(
                "a crossing at a singular point of a spline's curve",
            ))?;
            if sigma.cmp(&nil) == Ordering::Less || sigma.cmp(&one) == Ordering::Greater {
                continue;
            }
            if at[0].cmp(&peval(&qa.x, &sigma)) != Ordering::Equal
                || at[1].cmp(&peval(&qa.y, &sigma)) != Ordering::Equal
            {
                return Err(Error::ComputationLimit(
                    "a crossing off its spline's parametrization",
                ));
            }
            let (at0, at1) = (
                sigma.cmp(&nil) == Ordering::Equal,
                sigma.cmp(&one) == Ordering::Equal,
            );
            if at0 && k > 0 {
                // The previous arc's end: found there.
                continue;
            }
            if root.mult > 1 {
                return Err(tangent_splines());
            }
            let q_knot = at1 && k + 1 < n;
            let q_end = (at0 && k == 0) || (at1 && k + 1 == n);
            let tau_q = q.tau_of(k, &sigma.a);
            if root.knot && q_knot {
                return Err(touching_knot());
            }
            if root.knot {
                // q's tangent there, in p's frame.
                let dq = q.deriv(&tau_q);
                let t = apply_dir(&back, &dq);
                let (arrive, leave) = knot_legs(p, &root.tau);
                if !crosses_at_knot(&arrive, &leave, &t) {
                    return Err(touching_knot());
                }
            }
            if q_knot {
                let dp = p.deriv(&root.tau);
                let t = apply_dir(&m, &dp);
                let (arrive, leave) = knot_legs(q, &tau_q);
                if !crosses_at_knot(&arrive, &leave, &t) {
                    return Err(touching_knot());
                }
            }
            out.push(Crossing {
                tau: [root.tau, tau_q],
                uv,
                end: [root.end, q_end],
            });
        }
    }
    Ok(out)
}

/// A spline wall's section by another prism's spline wall: generatrices
/// over the crossings strictly inside both runs.
pub(super) fn spline_spline(
    px: &Prism,
    sx: &Arc<SplineSeg>,
    py: &Prism,
    sy: &Arc<SplineSeg>,
) -> Result<Section> {
    let cs = crossings(sx, &px.f, sy, &py.f)?;
    Ok(Section::Curves(
        cs.into_iter()
            .filter(|c| !c.end[0] && !c.end[1])
            .map(|c| Crv::Line {
                p: qpoint(&px.f, &c.uv[0], &c.uv[1], &Qd::rat(zero())),
                d: px.f.n.clone(),
            })
            .collect(),
    ))
}

/// Where a curve over a spline segment (a cap edge or a crease) meets
/// another prism's spline wall (segment `seg` on frame `f`).
pub(super) fn wallcrv_wall(curve: &WallCrv, f: &Affine, seg: &SplineSeg) -> Result<EdgeMeet> {
    let cs = crossings(&curve.seg, &curve.f, seg, f)?;
    Ok(EdgeMeet::Points(
        cs.into_iter()
            .map(|c| {
                let x = curve.point(&c.tau[0]);
                (Pos::T(c.tau[0].clone()), x)
            })
            .collect(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::SplineSpan;
    use crate::{BSplineCurve2, Frame3, Point2, Point3, Tolerance, Vec3};
    use num_bigint::BigInt;

    fn r(a: i64, b: i64) -> R {
        R::new(BigInt::from(a), BigInt::from(b))
    }

    fn segment(poles: &[(f64, f64)], op: usize) -> SplineSeg {
        let n = poles.len();
        let curve = BSplineCurve2::new(
            n - 1,
            poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect(),
            None,
            vec![0.0, 1.0],
            vec![n, n],
        )
        .expect("a Bezier arc");
        SplineSeg::new(&SplineSpan::whole(curve), false, op).expect("a segment")
    }

    fn affine(origin: (f64, f64, f64), normal: (f64, f64, f64), x: (f64, f64, f64)) -> Affine {
        let f = Frame3::new(
            Point3::new(origin.0, origin.1, origin.2),
            Vec3::new(normal.0, normal.1, normal.2),
            Vec3::new(x.0, x.1, x.2),
            Tolerance::default(),
        )
        .expect("a frame");
        Affine::new(&f).expect("its map")
    }

    /// The dome `y = x (4 - x) / 2`: its implicit equation vanishes on it
    /// (and not beside it), and its inversion gives back the parameter; a
    /// cubic of degree three in both coordinates likewise (the first
    /// subresultant's case).
    #[test]
    fn an_arc_s_implicit_equation_and_inversion() {
        let seg = segment(&[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)], 0);
        let imp = seg.arcs[0].implicit();
        for s in [r(0, 1), r(1, 3), r(1, 2), r(7, 9), r(1, 1)] {
            let p = seg.point(&Qd::rat(s.clone()));
            assert_eq!(imp.value(&p).sign(), Ordering::Equal);
            let back = imp.param(&p).expect("a parameter");
            assert_eq!(back.cmp(&Qd::rat(s)), Ordering::Equal);
        }
        let off = [Qd::rat(r(2, 1)), Qd::rat(r(1, 1))];
        assert_ne!(imp.value(&off).sign(), Ordering::Equal);
        let seg = segment(&[(0.0, 0.0), (1.0, 3.0), (3.0, -1.0), (4.0, 2.0)], 0);
        let imp = seg.arcs[0].implicit();
        for s in [r(1, 5), r(1, 2), r(5, 6)] {
            let p = seg.point(&Qd::rat(s.clone()));
            assert_eq!(imp.value(&p).sign(), Ordering::Equal);
            assert_eq!(
                imp.param(&p).expect("a parameter").cmp(&Qd::rat(s)),
                Ordering::Equal
            );
        }
    }

    /// The dome and the dome turned a quarter turn about the axis cross
    /// once, in the object's field, the same point found from either side.
    #[test]
    fn two_arcs_cross_in_the_object_s_field_either_way() {
        let a = segment(&[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)], 0);
        let b = segment(&[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)], 1);
        let fa = affine((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let fb = affine((3.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0));
        let ab = crossings(&a, &fa, &b, &fb).expect("crossings");
        let ba = crossings(&b, &fb, &a, &fa).expect("crossings");
        assert_eq!(ab.len(), 1);
        assert_eq!(ba.len(), 1);
        for k in 0..2 {
            assert_eq!(ab[0].tau[k].cmp(&ba[0].tau[1 - k]), Ordering::Equal);
        }
        // Both in the object's arc's field: one generator.
        let (g, h) = (ab[0].tau[0].gen(), ba[0].tau[1].gen());
        assert!(g.zip(h).is_some_and(|(g, h)| Arc::ptr_eq(g, h)));
        // The point lies on both curves, found again from its coordinates.
        let uv = &ab[0].uv;
        let tau_b = b
            .locate(&apply(&map2(&fa, &fb).expect("parallel"), uv))
            .expect("on the other");
        assert_eq!(tau_b.cmp(&ab[0].tau[1]), Ordering::Equal);
        let tau_a = a.locate(uv).expect("on its own");
        assert_eq!(tau_a.cmp(&ab[0].tau[0]), Ordering::Equal);
    }

    /// The dome against itself mirrored (its axis reversed, `x` along
    /// `-x`): one curve in both frames, refused as coincident walls; a
    /// crossing axis refused.
    #[test]
    fn coincident_arcs_and_crossing_axes_are_refused() {
        let a = segment(&[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)], 0);
        let b = segment(&[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)], 1);
        let fa = affine((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let mirror = affine((4.0, 0.0, 6.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0));
        assert!(matches!(
            crossings(&a, &fa, &b, &mirror),
            Err(Error::OutOfDomain(m)) if m.contains("one surface")
        ));
        let side = affine((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
        assert!(matches!(
            crossings(&a, &fa, &b, &side),
            Err(Error::OutOfDomain(m)) if m.contains("crossing axes")
        ));
    }
}
