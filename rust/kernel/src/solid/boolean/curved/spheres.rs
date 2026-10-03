//! S9d.2: spheres against prisms with arcs, and two spheres (REVIEW_NOTES.md,
//! S9d.2 refined). Two spheres meet in their radical plane's section. A
//! cylinder and a sphere meet where the cylinder's rulings meet the sphere:
//! over the cylinder's angle, the ruling's quadratic has a discriminant
//! positive all round (two rings over the angle, graphs of `Curve3::Meet`
//! with the sphere as the other quadric), negative all round (apart), or
//! changing sign (a loop, whose turning points need graphs over the height:
//! S9d.2b's, and in a turned frame S9d.2c's, `spheres_turned.rs`). A
//! prism's cap circle meets a sphere at the roots of a quartic in its
//! half-angle tangent (`algebraic.rs`); a sphere's own circle (a rim, the
//! split) meets a cylinder where `F0 + s F1 = 0`, `s` its radius over its
//! basis's length, at algebraic points with one surd (on a basis of
//! unequal axes, a turned cap's, S9d.2c's resultant), and another sphere on
//! the spheres' radical plane.
use super::meet::{tangency, CylPair, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::procedural::{other_of, other_sphere, MeetCrv, Other, Quartic, Ruled};
use super::sphere::{plane_section, Circ};
use super::turned::{middle, roots, square_sum, trim, Chart, Lin};
use crate::polynomial::real::IntPolynomial;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// How a cylinder (of operand `k`) and a sphere meet: rings over the
/// cylinder's angle, or apart.
pub(super) fn sphere_cyl(k: usize, cyl: Cyl, c: &V, r: &R, res: f64) -> Result<CylPair> {
    let other = other_sphere(c, r);
    let (a, d) = super::turned::discriminant(cyl, &other);
    // A discriminant vanishing identically: the cylinder tangent to the
    // sphere all round a circle (a coaxial pair of equal radii).
    let unit = Chart {
        c0: int(1),
        s0: zero(),
    };
    if trim(d.poly(&unit)).is_empty() && d.value(&[int(-1), zero()]) == zero() {
        return Err(tangency());
    }
    let chart = super::turned::negative_chart(&d)?;
    let test = chart.clone().unwrap_or(Chart {
        c0: int(1),
        s0: zero(),
    });
    super::turned::near_node(&a, &d, &test, res)?;
    let Some(chart) = chart else {
        // D >= 0 all round (and no double root): two rings.
        let piece = |plus| MeetCrv::new(k, cyl, &other, plus, None);
        return Ok(CylPair::Quartic(Box::new(Quartic {
            pieces: vec![piece(true), piece(false)],
            switches: Vec::new(),
        })));
    };
    if roots(&d.poly(&chart))?.is_empty() {
        return Ok(CylPair::Apart);
    }
    let slope = zero();
    let ruled = Ruled {
        f: cyl.0,
        c: cyl.1,
        r: cyl.2,
        k: &slope,
    };
    loops(k, ruled, c, r, &other, &d, &chart)
}

/// A piece of a cylinder's or a cone's meeting with a sphere over the
/// height (S9d.2b, S9d.3b.2): at height `w` the carrier's circle, of radius
/// `rho(w) = r + k w`, meets the sphere where `alpha cos u + beta sin u =
/// g(w) / rho(w)` (`alpha`, `beta` per unit radius), on the branch `plus`
/// (the sign of `alpha sin u - beta cos u`), over `range` (ascending
/// heights). On a turned cylinder (S9d.2c, `exact` false) the circle is an
/// ellipse in the world and the piece is the one root of the sphere's
/// function on it on that branch (`spheres_turned.rs`). Heights, branches
/// and the carrier's membership are read in the carrier's local coordinates
/// (`rows`, the inverse frame's: its axes in an exact frame).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct RiseCrv {
    pub(super) carrier: usize,
    pub(super) o: V,
    pub(super) x: V,
    pub(super) y: V,
    pub(super) n: V,
    /// The carrier frame's inverse rows (local `u`, `v`, `w`).
    rows: [V; 3],
    /// Whether the carrier's stored axes are exactly orthonormal.
    exact: bool,
    pub(super) r: R,
    /// The carrier's slope (zero for a cylinder).
    pub(super) k: R,
    pub(super) c: V,
    pub(super) rr: R,
    pub(super) plus: bool,
    pub(super) range: [Qd; 2],
}

impl RiseCrv {
    /// `alpha`, `beta` (per unit radius) and `g`'s coefficients `[g0, g1,
    /// g2]` (`rho(w)^2` expanded in them; exact in an exact frame, the
    /// half-turns' axis and a binary64 guide in a turned one).
    pub(super) fn coefficients(&self) -> (R, R, [R; 3]) {
        let d = sub(&self.o, &self.c);
        let (r, k) = (&self.r, &self.k);
        (
            int(2) * dot(&self.x, &d),
            int(2) * dot(&self.y, &d),
            [
                &self.rr * &self.rr - dot(&d, &d) - r * r,
                int(-2) * dot(&self.n, &d) - int(2) * r * k,
                int(-1) - k * k,
            ],
        )
    }

    /// The carrier's ruled surface and its quadric, and the sphere (S9e.3b:
    /// the curve's two surfaces).
    pub(super) fn surfaces(&self) -> (super::triple::Ruling, Other, Other) {
        let ruling = super::triple::Ruling {
            o: self.o.clone(),
            x: self.x.clone(),
            y: self.y.clone(),
            n: self.n.clone(),
            r: self.r.clone(),
            k: self.k.clone(),
        };
        let carrier = Other {
            g: vec![self.rows[0].clone(), self.rows[1].clone()],
            e: vec![dot(&self.rows[0], &self.o), dot(&self.rows[1], &self.o)],
            r: self.r.clone(),
            h: self.rows[2].clone(),
            eh: dot(&self.rows[2], &self.o),
            t: self.k.clone(),
        };
        (ruling, carrier, other_sphere(&self.c, &self.rr))
    }

    /// The carrier's radius at a height.
    fn rho(&self, w: &R) -> R {
        &self.r + &self.k * w
    }

    /// The point at a rational height (`None` off the piece's branch
    /// domain).
    pub(super) fn at(&self, w: &R) -> Result<Option<QV>> {
        if !self.exact {
            return super::spheres_turned::rise_at(self, w);
        }
        let (a, b, g) = self.coefficients();
        let rho = self.rho(w);
        if rho <= zero() {
            return Ok(None);
        }
        let gw = (&g[0] + &g[1] * w + &g[2] * w * w) / &rho;
        let Some(sols) = super::meet::trig(&a, &b, &gw)? else {
            return Ok(None);
        };
        let Some(cs) = sols.get(usize::from(self.plus)) else {
            return Ok(None);
        };
        let p = qadd(
            &qadd(&qv(&self.o), &qscale(&self.x, &cs[0].scale(&rho))),
            &qadd(
                &qscale(&self.y, &cs[1].scale(&rho)),
                &qv(&scale(&self.n, w)),
            ),
        );
        Ok(Some(p))
    }

    /// The carrier's circle centre at height 0.
    pub(super) fn o_model(&self) -> V {
        self.o.clone()
    }

    /// A point's local coordinates from the circle's centre at height 0.
    fn local(&self, p: &QV) -> QV {
        let d = qsub(p, &qv(&self.o));
        [0, 1, 2].map(|k| qdot(&d, &self.rows[k]))
    }

    /// A point's height.
    pub(super) fn height(&self, p: &QV) -> Qd {
        qdot(&qsub(p, &qv(&self.o)), &self.rows[2])
    }

    fn branch(&self, p: &QV) -> Ordering {
        let (a, b, _) = self.coefficients();
        let [lu, lv, _] = self.local(p);
        lv.scale(&a).sub(&lu.scale(&b)).sign()
    }

    /// Whether a point lies on the piece (both surfaces, its branch, its
    /// heights, ends included).
    pub(super) fn on(&self, p: &QV) -> bool {
        let [dx, dy, dw] = self.local(p);
        let rho = dw.scale(&self.k).add_r(&self.r);
        let cyl = dx.mul(&dx).add(&dy.mul(&dy)).sub(&rho.mul(&rho));
        let e = qsub(p, &qv(&self.c));
        let sph = qqdot(&e, &e).add_r(&-(&self.rr * &self.rr));
        if cyl.sign() != Ordering::Equal || sph.sign() != Ordering::Equal {
            return false;
        }
        let want = if self.plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if self.branch(p) != want {
            return false;
        }
        let h = self.height(p);
        h.cmp(&self.range[0]) != Ordering::Less && h.cmp(&self.range[1]) != Ordering::Greater
    }

    /// The unit-free tangent at a point, running up.
    pub(super) fn tangent(&self, p: &QV) -> QV {
        let [du, dv, dw] = self.local(p);
        // The carrier's normal (a cone's leans against its axis): the
        // gradient of `u^2 + v^2 - rho(w)^2` in the world, halved.
        let rho = dw.scale(&self.k).add_r(&self.r);
        let gc = qadd(
            &qadd(&qscale(&self.rows[0], &du), &qscale(&self.rows[1], &dv)),
            &qscale(&self.rows[2], &rho.scale(&-self.k.clone())),
        );
        let t = qcross(&gc, &qsub(p, &qv(&self.c)));
        if qdot(&t, &self.rows[2]).sign() == Ordering::Less {
            t.map(|x| x.neg())
        } else {
            t
        }
    }

    /// Binary64 points from height `w0` to `w1`.
    pub(super) fn samples(&self, w0: f64, w1: f64, n: usize) -> Vec<[f64; 3]> {
        let f = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (o, x, y, nn) = (f(&self.o), f(&self.x), f(&self.y), f(&self.n));
        let (a, b, g) = self.coefficients();
        let (a, b) = (rational_f64(&a), rational_f64(&b));
        let g = g.map(|x| rational_f64(&x));
        let (r, k) = (rational_f64(&self.r), rational_f64(&self.k));
        let sign = if self.plus { 1.0 } else { -1.0 };
        (0..=n)
            .map(|i| {
                let w = w0 + (w1 - w0) * i as f64 / n as f64;
                let rho = r + k * w;
                let q = ((g[0] + g[1] * w + g[2] * w * w) / (rho * a.hypot(b))).clamp(-1.0, 1.0);
                let u = b.atan2(a) + sign * q.acos();
                [0, 1, 2].map(|j| o[j] + rho * (u.cos() * x[j] + u.sin() * y[j]) + w * nn[j])
            })
            .collect()
    }
}

/// A sphere and a cylinder meeting in loops (S9d.2b), the cylinder's frame
/// exact: each interval of the cylinder's discriminant's positive values
/// holds a loop, its turning points of both kinds (the ruling tangent, the
/// circle tangent) ordered along it and rational switches between those of
/// different kinds: graphs over the height about the first, over the angle
/// about the second, each verified exactly. On a turned cylinder (S9d.2c)
/// or cone (S9d.3c) the height graph's turning points are the
/// discriminant's real roots (`spheres_turned::Height`) and its pieces are
/// verified by its checks.
pub(super) fn loops(
    k: usize,
    ruled: Ruled,
    c: &V,
    rr: &R,
    other: &Other,
    d: &super::turned::Form,
    chart: &Chart,
) -> Result<CylPair> {
    let Ruled {
        f,
        c: cc,
        r,
        k: slope,
    } = ruled;
    let exact = f.orthonormal();
    let o = f.point(&cc[0], &cc[1], &zero());
    let rise = |plus: bool, range: [Qd; 2]| RiseCrv {
        carrier: k,
        o: o.clone(),
        x: f.x.clone(),
        y: f.y.clone(),
        n: f.n.clone(),
        rows: [f.row(0).clone(), f.row(1).clone(), f.row(2).clone()],
        exact,
        r: r.clone(),
        k: slope.clone(),
        c: c.clone(),
        rr: rr.clone(),
        plus,
        range,
    };
    let probe = rise(true, [Qd::rat(zero()), Qd::rat(zero())]);
    let (alpha, beta, g) = probe.coefficients();
    let rho2 = &alpha * &alpha + &beta * &beta;
    if rho2 == zero() {
        return Err(Error::ComputationLimit(
            "a loop about a sphere centred on the carrier's axis",
        ));
    }
    // S9d.2c: on a turned cylinder (S9d.3c: or cone) the exact height graph
    // and its checks.
    let height = if exact {
        None
    } else {
        Some(super::spheres_turned::Height::new(
            (&o, &f.x, &f.y, &f.n),
            r,
            slope,
            c,
            rr,
            [alpha.clone(), beta.clone()],
        )?)
    };
    // The height graph's discriminant rho^2 rho(w)^2 - g(w)^2, a quartic in
    // w.
    let gp = trim(g.to_vec());
    let rw = trim(vec![r.clone(), slope.clone()]);
    let dw = trim(padd(
        &pscale(&pmul(&rw, &rw), &rho2),
        &pscale(&pmul(&gp, &gp), &int(-1)),
    ));
    let meet =
        |plus: bool, range: Option<[[Qd; 2]; 2]>| MeetCrv::ruled(k, ruled, other, plus, range);
    let pa = d.poly(chart);
    let mut ra = roots(&pa)?;
    let mut rws = match &height {
        Some(h) => h.turns(),
        None => roots(&dw)?,
    };
    for r in ra.iter_mut().chain(rws.iter_mut()) {
        r.refine_for_signs(160);
    }
    if ra.len() % 2 != 0 {
        return Err(Error::ComputationLimit(
            "an odd count of a cylinder's turning points",
        ));
    }
    // Turning points of the height graph (g = +-rho: u = phi or phi + pi):
    // their chart t and the angle graph's branch there, binary64 views.
    let (af, bf) = (rational_f64(&alpha), rational_f64(&beta));
    let phi = bf.atan2(af);
    let fl = |v: &V| v.clone().map(|y| rational_f64(&y));
    let (of, xf, yf, nf, cf) = (fl(&o), fl(&f.x), fl(&f.y), fl(&f.n), fl(c));
    let (rf, kf) = (rational_f64(r), rational_f64(slope));
    let mut w_turns: Vec<(f64, bool)> = Vec::new();
    for root in &rws {
        let w = rational_f64(&middle(root));
        let rf = rf + kf * w;
        let gw = rational_f64(&g[0]) + rational_f64(&g[1]) * w + rational_f64(&g[2]) * w * w;
        // The carrier's angle where `alpha cos + beta sin = g / rho`: on a
        // cone's far nappe (`rho < 0`, S9d.3c) the other side of `phi`.
        let u = if gw * rf > 0.0 {
            phi
        } else {
            phi + std::f64::consts::PI
        };
        let p: [f64; 3] =
            [0, 1, 2].map(|j| of[j] + rf * (u.cos() * xf[j] + u.sin() * yf[j]) + w * nf[j]);
        // The branch: the sphere's gradient along the ruling there (a cone's
        // leaning out, S9d.3c; the axis for a cylinder).
        let dir: [f64; 3] = [0, 1, 2].map(|j| nf[j] + kf * (u.cos() * xf[j] + u.sin() * yf[j]));
        let up = (0..3).map(|j| (p[j] - cf[j]) * dir[j]).sum::<f64>();
        if up.abs() < 1e-12 * (1.0 + rf.abs()) {
            return Err(Error::ComputationLimit("a turning point of both graphs"));
        }
        // The chart's t of (cos u, sin u): rotated back by the base.
        let (c0, s0) = (rational_f64(&chart.c0), rational_f64(&chart.s0));
        let (cu, su) = (u.cos(), u.sin());
        let (cr, sr) = (c0 * cu + s0 * su, c0 * su - s0 * cu);
        w_turns.push((sr / (1.0 + cr), up > 0.0));
    }
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for pair in ra.chunks(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let (a64, b64) = (rational_f64(&middle(a)), rational_f64(&middle(b)));
        // Events along the loop: lambda in [0, 2), the + branch (the angle
        // graph's) from a to b, then the - branch back; true for the
        // ruling's turning points.
        let mut events: Vec<(f64, bool)> = vec![(0.0, true), (1.0, true)];
        for &(t, plus) in &w_turns {
            if t <= a64 || t >= b64 {
                continue;
            }
            let l = (t - a64) / (b64 - a64);
            events.push((if plus { l } else { 2.0 - l }, false));
        }
        events.sort_by(|x, y| x.0.total_cmp(&y.0));
        let n = events.len();
        if n < 3 {
            return Err(Error::ComputationLimit(
                "a loop without turning points of the other graph",
            ));
        }
        for i in 0..n {
            let l1 = events[(i + 1) % n].0 + if i + 1 == n { 2.0 } else { 0.0 };
            if l1 - events[i].0 < 1e-9 {
                return Err(Error::Degenerate("two turning points within rounding"));
            }
        }
        let at = |l: f64| -> (R, bool) {
            let l = l.rem_euclid(2.0);
            if l < 1.0 {
                (q(a64 + l * (b64 - a64)), true)
            } else {
                (q(b64 - (l - 1.0) * (b64 - a64)), false)
            }
        };
        let point = |l: f64| -> Result<(QV, R, bool)> {
            let (t, plus) = at(l);
            if a.compare_rational(&t) != Ordering::Less
                || b.compare_rational(&t) != Ordering::Greater
            {
                return Err(Error::ComputationLimit("a switch point off its loop"));
            }
            let p = meet(plus, None)
                .at(&chart.at(&t))
                .ok_or(Error::ComputationLimit("a switch point off its piece"))?;
            Ok((p, t, plus))
        };
        let mut sw: Vec<(f64, QV, R, bool)> = Vec::new();
        for i in 0..n {
            let (e0, e1) = (events[i], events[(i + 1) % n]);
            if e0.1 == e1.1 {
                continue;
            }
            let l1 = e1.0 + if i + 1 == n { 2.0 } else { 0.0 };
            let l = 0.5 * (e0.0 + l1);
            let (p, t, plus) = point(l)?;
            sw.push((l.rem_euclid(2.0), p, t, plus));
        }
        sw.sort_by(|x, y| x.0.total_cmp(&y.0));
        let m = sw.len();
        if m < 2 {
            return Err(Error::ComputationLimit("a loop with one switch"));
        }
        for i in 0..m {
            let (s0, s1) = (&sw[i], &sw[(i + 1) % m]);
            let l1 = s1.0 + if i + 1 == m { 2.0 } else { 0.0 };
            let first = events
                .iter()
                .find(|e| {
                    let l = if e.0 < s0.0 { e.0 + 2.0 } else { e.0 };
                    l > s0.0 && l < l1
                })
                .ok_or(Error::ComputationLimit("a run without turning points"))?;
            if first.1 {
                // The ruling's turning points: a graph over the height.
                let (h0, h1) = (probe.height(&s0.1), probe.height(&s1.1));
                let (lo, hi) = if h0.cmp(&h1) == Ordering::Less {
                    (h0, h1)
                } else {
                    (h1, h0)
                };
                let sign = probe.branch(&s0.1);
                if sign == Ordering::Equal || probe.branch(&s1.1) != sign {
                    return Err(Error::ComputationLimit(
                        "a run over the height changing branch",
                    ));
                }
                let clear = match &height {
                    Some(h) => h.verified(&lo, &hi, sign == Ordering::Greater)?,
                    None => {
                        let chain = super::turned::sturm(&dw);
                        super::turned::changes(&chain, &lo) == super::turned::changes(&chain, &hi)
                    }
                };
                if !clear {
                    return Err(Error::ComputationLimit("a turning point inside a piece"));
                }
                pieces.push(Crv::Rise(Box::new(rise(
                    sign == Ordering::Greater,
                    [lo, hi],
                ))));
            } else {
                // The circle's turning points: a graph over the angle.
                if s0.3 != s1.3 {
                    return Err(Error::ComputationLimit(
                        "a run over the angle changing branch",
                    ));
                }
                let (lo, hi) = if s0.2 < s1.2 {
                    (&s0.2, &s1.2)
                } else {
                    (&s1.2, &s0.2)
                };
                if ra.iter().any(|r| {
                    r.compare_rational(lo) == Ordering::Greater
                        && r.compare_rational(hi) == Ordering::Less
                }) {
                    return Err(Error::ComputationLimit("a turning point inside a piece"));
                }
                let rng = [lo, hi].map(|t| {
                    let cs = chart.at(t);
                    [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]
                });
                pieces.push(Crv::Meet(Box::new(meet(s0.3, Some(rng)))));
            }
        }
        switches.extend(sw.into_iter().map(|s| s.1));
    }
    Ok(CylPair::Mixed(Box::new(Mixed { pieces, switches })))
}

/// A cylinder's and a sphere's loops: pieces over the angle and over the
/// height, and their switches (S9d.2b).
#[derive(Debug, Clone)]
pub(super) struct Mixed {
    pub(super) pieces: Vec<Crv>,
    pub(super) switches: Vec<QV>,
}

/// Two spheres' meeting: their radical plane's section (`None` apart), the
/// same sphere refused, the plane crossing either sphere within the
/// resolution `res` of tangency `Degenerate`.
pub(super) fn sphere_sphere(c1: &V, r1: &R, c2: &V, r2: &R, res: f64) -> Result<Option<Circ>> {
    let m = sub(c2, c1);
    if is_zero(&m) {
        return Err(Error::Degenerate("two spheres about one centre"));
    }
    let (p0, m) = radical(c1, r1, c2, r2);
    plane_section(c2, r2, &p0, &m, res)?;
    plane_section(c1, r1, &p0, &m, res)
}

/// The radical plane of two spheres: a point and its normal.
pub(super) fn radical(c1: &V, r1: &R, c2: &V, r2: &R) -> (V, V) {
    let m = sub(c2, c1);
    let rhs = dot(c2, c2) - dot(c1, c1) + r1 * r1 - r2 * r2;
    let p0 = scale(&m, &(rhs / (int(2) * dot(&m, &m))));
    (p0, m)
}

/// Where a sphere's own circle meets a cylinder (`f`, `cy`, `ry`): points
/// `c + s (cos x + sin y)`, `s^2 = r2 / |x|^2` (a basis of equal lengths),
/// where `F0 + s F1 = 0`; S9d.2c's resultant for another basis (`res` the
/// resolution its tangencies are refused within).
pub(super) fn circ_cylinder(
    circ: &Circ,
    f: &Affine,
    cy: &P2,
    ry: &R,
    res: f64,
) -> Result<EdgeMeet> {
    circ_quadric(circ, &other_of(f, cy, ry), res)
}

/// Whether a sphere's circle lies on a basis of equal axes (an exact frame's
/// or a whole sphere's; a turned cap's are unequal, S9d.2c).
pub(super) fn equal_axes(circ: &Circ) -> bool {
    dot(&circ.x, &circ.x) == dot(&circ.y, &circ.y) && dot(&circ.x, &circ.y) == zero()
}

pub(super) fn circ_quadric(circ: &Circ, o: &Other, res: f64) -> Result<EdgeMeet> {
    if !equal_axes(circ) {
        return super::spheres_turned::circ_ellipse(circ, o, res);
    }
    let xx = dot(&circ.x, &circ.x);
    let sigma = &circ.r2 / &xx;
    // a_i + s L_i(cos, sin): F0 = sum a_i^2 + sigma sum L_i^2 - R^2, F1 = 2
    // sum a_i L_i.
    let a: Vec<R> = (0..o.g.len())
        .map(|i| dot(&o.g[i], &circ.c) - &o.e[i])
        .collect();
    let l: Vec<Lin> = (0..o.g.len())
        .map(|i| [zero(), dot(&o.g[i], &circ.x), dot(&o.g[i], &circ.y)])
        .collect();
    // The other's radius term `ra + s rl` (a cone's, S9d.3b; `r` alone
    // otherwise) squared and subtracted.
    let rad = o.radius_lin(&circ.c, &circ.x, &circ.y);
    let (ra, rl): (R, Lin) = (rad[0].clone(), [zero(), rad[1].clone(), rad[2].clone()]);
    let mut f0 = square_sum(&l)
        .sub(&square_sum(std::slice::from_ref(&rl)))
        .scaled(&sigma);
    f0.add_const(&(a.iter().fold(zero(), |acc, x| acc + x * x) - &ra * &ra));
    let f1: Lin = [0, 1, 2].map(|j| {
        a.iter()
            .zip(&l)
            .fold(zero(), |acc, (ai, li)| acc + int(2) * ai * &li[j])
            - int(2) * &ra * &rl[j]
    });
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let s = rational_sqrt(&sigma);
    // F0 times (1 + t^2)^2, F1 times (1 + t^2).
    let p0 = trim(f0.poly(&chart));
    let [cn, sn] = chart.numerators();
    let w = vec![int(1), zero(), int(1)];
    let p1 = trim(padd(
        &padd(&pscale(&w, &f1[0]), &pscale(&cn, &f1[1])),
        &pscale(&sn, &f1[2]),
    ));
    // The equation's polynomial: F0 + s F1 (s rational), or F0^2 - sigma
    // F1^2 (its roots of the right signs kept).
    let eq = match &s {
        Some(s) => trim(padd(&p0, &pmul(&pscale(&p1, s), &w))),
        None => trim(padd(
            &pmul(&p0, &p0),
            &pscale(&pmul(&pmul(&p1, &p1), &pmul(&w, &w)), &-sigma.clone()),
        )),
    };
    if eq.is_empty() {
        return Err(tangency());
    }
    let point = |cs: &[K; 2]| -> QV {
        let dir = [0, 1, 2].map(|j| cs[0].scale(&circ.x[j]).add(&cs[1].scale(&circ.y[j])));
        [0, 1, 2].map(|j| match &s {
            Some(s) => Qd::of(K::Rat(circ.c[j].clone()).add(&dir[j].scale(s))),
            None => Qd::parts(K::Rat(circ.c[j].clone()), dir[j].clone(), sigma.clone()),
        })
    };
    let place = |cs: &[K; 2]| -> [Qd; 2] {
        match &s {
            Some(s) => [Qd::of(cs[0].scale(s)), Qd::of(cs[1].scale(s))],
            None => [
                Qd::parts(K::Rat(zero()), cs[0].clone(), sigma.clone()),
                Qd::parts(K::Rat(zero()), cs[1].clone(), sigma.clone()),
            ],
        }
    };
    let mut out = Vec::new();
    // The chart's antipode (-1, 0): where the equation drops degrees.
    let anti = [K::Rat(int(-1)), K::Rat(zero())];
    let f0a = f0.value(&[int(-1), zero()]);
    let f1a = &f1[0] - &f1[1];
    let at_anti = match &s {
        Some(s) => f0a.clone() + s * &f1a == zero(),
        None => &f0a * &f0a == &sigma * &f1a * &f1a && (f0a == zero() || sign(&f0a) != sign(&f1a)),
    };
    if at_anti {
        out.push((Pos::Ang(place(&anti)), point(&anti)));
    }
    let deg_ok = |r: &crate::polynomial::real::AlgebraicRoot| -> bool {
        match &s {
            Some(_) => true,
            None => {
                let (g0, g1) = (
                    r.sign_polynomial(&IntPolynomial::from_rationals(&p0)),
                    r.sign_polynomial(&IntPolynomial::from_rationals(&p1)),
                );
                g0 == Ordering::Equal || (g1 != Ordering::Equal && g0 != g1)
            }
        }
    };
    for root in roots(&eq)? {
        if !deg_ok(&root) {
            continue;
        }
        let g = Arc::new(Gen::new(eq.clone(), root));
        let t = K::generator(&g);
        let den = t.mul(&t).add(&K::Rat(int(1)));
        let inv = den
            .recip()
            .ok_or(Error::ComputationLimit("a circle's crossing at infinity"))?;
        let cs = [
            K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv),
            t.scale(&int(2)).mul(&inv),
        ];
        out.push((Pos::Ang(place(&cs)), point(&cs)));
    }
    Ok(EdgeMeet::Points(out))
}

fn padd(a: &[R], b: &[R]) -> Vec<R> {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| {
                a.get(i).cloned().unwrap_or_else(zero) + b.get(i).cloned().unwrap_or_else(zero)
            })
            .collect(),
    )
}

fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

fn pscale(a: &[R], k: &R) -> Vec<R> {
    trim(a.iter().map(|x| x * k).collect())
}

/// Where a sphere's own circle (on the sphere `(c1, r1)`) meets another
/// sphere: on the two spheres' radical plane.
pub(super) fn circ_sphere(circ: &Circ, c1: &V, r1: &R, c2: &V, r2: &R) -> Result<EdgeMeet> {
    if is_zero(&sub(c2, c1)) {
        return Err(Error::Degenerate("two spheres about one centre"));
    }
    let (p0, m) = radical(c1, r1, c2, r2);
    Ok(match circ.meet_plane(&p0, &m)? {
        None => EdgeMeet::Along,
        Some(cs) => EdgeMeet::Points(
            cs.into_iter()
                .map(|e| {
                    let x = qadd(
                        &qv(&circ.c),
                        &qadd(&qscale(&circ.x, &e[0]), &qscale(&circ.y, &e[1])),
                    );
                    (Pos::Ang(e), x)
                })
                .collect(),
        ),
    })
}
