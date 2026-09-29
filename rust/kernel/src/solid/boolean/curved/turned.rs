//! S9c.2b.1: two cylinders in turned frames (affine models on their stored
//! axes, elliptic in the world) whose axes cross, meeting in a quartic
//! (REVIEW_NOTES.md, S9c.2b refined). Each cylinder's discriminant `D_K`
//! over its angle is a quartic in the half-angle tangent `t` of a chart,
//! its roots (the turning points of `K`'s graph) isolated exactly. `D_K > 0`
//! all round gives two rings over `K`; otherwise every interval of
//! `D_A >= 0` holds a loop, its turning points of both kinds ordered along
//! it and a rational point of `A`'s angle placed between each adjacent pair
//! of different kinds; the pieces are then verified exactly (no root of the
//! carrier's discriminant in a piece's range). A double root, or an
//! extremum of `D_K` within the resolution of zero (two branches that
//! close), is `Degenerate`. A cap's circle meets the other cylinder where a
//! quartic in its own `t` vanishes: roots on the edge's arc are S9c.2b.2's.
use super::graph::between_ccw;
use super::meet::{tangency, CylPair};
use super::model::*;
use super::num::*;
use super::procedural::{other_of, MeetCrv, Other, Quartic};
use crate::polynomial::real::{isolate, AlgebraicRoot, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A section crossing a cap's circle in turned frames: an algebraic point.
pub(super) fn algebraic() -> Error {
    Error::OutOfDomain("two cylinders' section crossing a cap's circle in turned frames (S9c.2b.2)")
}

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

// ------------------------------------------------------------ polynomials

/// A polynomial with rational coefficients, ascending powers.
type Poly = Vec<R>;

fn trim(mut p: Poly) -> Poly {
    while p.last().is_some_and(|c| *c == zero()) {
        p.pop();
    }
    p
}

fn padd(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| {
                a.get(i).cloned().unwrap_or_else(zero) + b.get(i).cloned().unwrap_or_else(zero)
            })
            .collect(),
    )
}

fn pmul(a: &Poly, b: &Poly) -> Poly {
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

fn pscale(a: &Poly, k: &R) -> Poly {
    trim(a.iter().map(|x| x * k).collect())
}

fn pderiv(a: &Poly) -> Poly {
    trim(
        a.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * int(i as i64))
            .collect(),
    )
}

/// The remainder of `a` by `b` (`b` nonzero).
fn prem(a: &Poly, b: &Poly) -> Poly {
    let mut r = a.clone();
    let lead = b.last().expect("a nonzero divisor").clone();
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let c = r.last().expect("nonempty") / &lead;
        for (i, x) in b.iter().enumerate() {
            r[i + shift] -= &c * x;
        }
        r = trim(r);
    }
    r
}

/// The Sturm chain of a square-free polynomial.
fn sturm(p: &Poly) -> Vec<Poly> {
    let mut chain = vec![p.clone(), pderiv(p)];
    while !chain.last().expect("a chain").is_empty() {
        let n = chain.len();
        let r = prem(&chain[n - 2], &chain[n - 1]);
        if r.is_empty() {
            break;
        }
        chain.push(pscale(&r, &int(-1)));
    }
    chain.retain(|x| !x.is_empty());
    chain
}

/// A polynomial's exact value at a surd.
fn eval(p: &Poly, x: &Qd) -> Qd {
    let mut acc = Qd::rat(zero());
    for c in p.iter().rev() {
        acc = acc.mul(x).add_r(c);
    }
    acc
}

/// Sign changes of a Sturm chain at a surd (zeros skipped).
fn changes(chain: &[Poly], x: &Qd) -> usize {
    let signs: Vec<Ordering> = chain
        .iter()
        .map(|p| eval(p, x).sign())
        .filter(|s| *s != Ordering::Equal)
        .collect();
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

fn int_poly(p: &Poly) -> IntPolynomial {
    IntPolynomial::from_rationals(p)
}

/// The real roots of a nonzero polynomial, exactly isolated in increasing
/// order; a repeated root is a tangency.
fn roots(p: &Poly) -> Result<Vec<AlgebraicRoot>> {
    let ip = int_poly(p);
    if ip.is_zero() {
        return Err(tangency());
    }
    if ip.is_constant() {
        return Ok(Vec::new());
    }
    if !ip.gcd(&ip.derivative()).is_constant() {
        return Err(tangency());
    }
    let abs = |x: &R| if *x < zero() { -x.clone() } else { x.clone() };
    let lead = abs(p.last().expect("nonzero"));
    let bound = p[..p.len() - 1]
        .iter()
        .map(|c| abs(c) / &lead)
        .fold(zero(), |m, x| if x > m { x } else { m })
        + int(1);
    isolate(
        &ip,
        -bound.clone(),
        bound,
        &mut Budget::new(RootIsolationOptions::default()),
    )
}

/// A root's midpoint, rational.
fn middle(r: &AlgebraicRoot) -> R {
    let (a, b) = r.isolator();
    (a + b) / int(2)
}

// ------------------------------------------------------------ charts

/// A chart of the circle: `(cos, sin)` the base `(c0, s0)` turned by
/// `2 atan t`, counter-clockwise with `t`; the base's antipode at infinity.
#[derive(Debug, Clone)]
struct Chart {
    c0: R,
    s0: R,
}

impl Chart {
    /// The numerators of `cos` and `sin` over `1 + t^2`.
    fn numerators(&self) -> [Poly; 2] {
        let (c0, s0) = (&self.c0, &self.s0);
        [
            trim(vec![c0.clone(), -(int(2) * s0), -c0.clone()]),
            trim(vec![s0.clone(), int(2) * c0, -s0.clone()]),
        ]
    }

    fn at(&self, t: &R) -> [R; 2] {
        let den = int(1) + t * t;
        let (c, s) = ((int(1) - t * t) / &den, int(2) * t / &den);
        [&self.c0 * &c - &self.s0 * &s, &self.s0 * &c + &self.c0 * &s]
    }

    /// The chart's `t` of a direction (`None` at the antipode).
    fn t_of(&self, cs: &[Qd; 2]) -> Option<Qd> {
        let c = cs[0].scale(&self.c0).add(&cs[1].scale(&self.s0));
        let s = cs[1].scale(&self.c0).sub(&cs[0].scale(&self.s0));
        let den = c.add_r(&int(1));
        if den.sign() == Ordering::Equal {
            return None;
        }
        // s / (1 + c) with den = a + b sqrt(d): times the conjugate.
        let conj = Qd::new(den.a.clone(), -den.b.clone(), den.d.clone());
        let norm = den.mul(&conj);
        debug_assert!(norm.is_rational());
        Some(s.mul(&conj).scale(&(int(1) / norm.a)))
    }
}

/// A quadratic form in `(cos, sin, 1)`: `cc c^2 + cs c s + ss s^2 + c1 c +
/// s1 s + k`.
#[derive(Debug, Clone, Default)]
struct Form {
    cc: R,
    cs: R,
    ss: R,
    c1: R,
    s1: R,
    k: R,
}

/// A linear form `l0 + l1 c + l2 s`.
type Lin = [R; 3];

fn square_sum(ls: &[Lin]) -> Form {
    let mut f = Form::default();
    for l in ls {
        f.cc += &l[1] * &l[1];
        f.cs += int(2) * &l[1] * &l[2];
        f.ss += &l[2] * &l[2];
        f.c1 += int(2) * &l[0] * &l[1];
        f.s1 += int(2) * &l[0] * &l[2];
        f.k += &l[0] * &l[0];
    }
    f
}

impl Form {
    fn scaled(&self, a: &R) -> Self {
        Self {
            cc: &self.cc * a,
            cs: &self.cs * a,
            ss: &self.ss * a,
            c1: &self.c1 * a,
            s1: &self.s1 * a,
            k: &self.k * a,
        }
    }

    fn sub(&self, o: &Self) -> Self {
        Self {
            cc: &self.cc - &o.cc,
            cs: &self.cs - &o.cs,
            ss: &self.ss - &o.ss,
            c1: &self.c1 - &o.c1,
            s1: &self.s1 - &o.s1,
            k: &self.k - &o.k,
        }
    }

    fn value(&self, cs: &[R; 2]) -> R {
        let (c, s) = (&cs[0], &cs[1]);
        &self.cc * c * c
            + &self.cs * c * s
            + &self.ss * s * s
            + &self.c1 * c
            + &self.s1 * s
            + &self.k
    }

    /// Times `(1 + t^2)^2` in a chart: a quartic in `t`.
    fn poly(&self, chart: &Chart) -> Poly {
        let [cn, sn] = chart.numerators();
        let w = vec![int(1), zero(), int(1)];
        let mut p = pscale(&pmul(&cn, &cn), &self.cc);
        p = padd(&p, &pscale(&pmul(&cn, &sn), &self.cs));
        p = padd(&p, &pscale(&pmul(&sn, &sn), &self.ss));
        p = padd(
            &p,
            &pmul(&padd(&pscale(&cn, &self.c1), &pscale(&sn, &self.s1)), &w),
        );
        padd(&p, &pscale(&pmul(&w, &w), &self.k))
    }
}

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// A ruling's quadratic against the other cylinder, `A w^2 + 2 B w + C`:
/// `A`, and the discriminant `B^2 - A C` as a form in the carrier's
/// `(cos, sin)`.
fn discriminant(k: Cyl, o: &Other) -> (R, Form) {
    let (f, c, r) = k;
    let base = f.point(&c[0], &c[1], &zero());
    let lin: Vec<Lin> = (0..2)
        .map(|i| {
            [
                dot(&o.g[i], &base) - &o.e[i],
                r * dot(&o.g[i], &f.x),
                r * dot(&o.g[i], &f.y),
            ]
        })
        .collect();
    let n: Vec<R> = (0..2).map(|i| dot(&o.g[i], &f.n)).collect();
    let a = &n[0] * &n[0] + &n[1] * &n[1];
    let b: Lin = [0, 1, 2].map(|j| &n[0] * &lin[0][j] + &n[1] * &lin[1][j]);
    let mut cform = square_sum(&lin);
    cform.k -= &o.r * &o.r;
    (a.clone(), square_sum(&[b]).sub(&cform.scaled(&a)))
}

/// A chart whose antipode has `D < 0` (a loop's piece never reaches it),
/// or `None` when `D >= 0` on the four axis points and every root gap.
fn negative_chart(d: &Form) -> Result<Option<Chart>> {
    for (c0, s0) in [(1, 0), (0, 1), (-1, 0), (0, -1)] {
        if d.value(&[int(-c0), int(-s0)]) < zero() {
            return Ok(Some(Chart {
                c0: int(c0),
                s0: int(s0),
            }));
        }
    }
    // A negative point between two roots of a chart with a nonzero lead.
    for (c0, s0) in [(1, 0), (0, 1)] {
        let chart = Chart {
            c0: int(c0),
            s0: int(s0),
        };
        if d.value(&[int(-c0), int(-s0)]) == zero() {
            continue;
        }
        let rs = roots(&d.poly(&chart))?;
        for w in rs.windows(2) {
            let t = (w[0].isolator().1 + w[1].isolator().0) / int(2);
            let p = chart.at(&t);
            if d.value(&p) < zero() {
                return Ok(Some(Chart {
                    c0: -p[0].clone(),
                    s0: -p[1].clone(),
                }));
            }
        }
        return Ok(None);
    }
    Err(limit(
        "a cylinders' discriminant vanishing at every axis point",
    ))
}

/// Two branches within the resolution: an extremum of `D` (a critical
/// point of its chart's quartic) where `|D| < (A res / 2)^2`, or a
/// repeated root (a tangency, `roots`).
fn near_node(a: &R, d: &Form, chart: &Chart, res: f64) -> Result<()> {
    let p = d.poly(chart);
    let _ = roots(&p)?;
    let delta = {
        let h = a * q(res) / int(2);
        &h * &h
    };
    let w2 = pmul(&vec![int(1), zero(), int(1)], &vec![int(1), zero(), int(1)]);
    let lo = int_poly(&padd(&p, &pscale(&w2, &delta)));
    let hi = int_poly(&padd(&p, &pscale(&w2, &-delta.clone())));
    let dp = pderiv(&p);
    if dp.is_empty() {
        return Ok(());
    }
    let crit = match roots(&dp) {
        Ok(r) => r,
        // A repeated critical point is an inflection, not a node.
        Err(Error::Degenerate(_)) => return Ok(()),
        Err(e) => return Err(e),
    };
    for r in crit {
        if r.sign_polynomial(&lo) == Ordering::Greater && r.sign_polynomial(&hi) == Ordering::Less {
            return Err(Error::Degenerate(
                "two cylinders' section within the resolution of a node",
            ));
        }
    }
    Ok(())
}

// ------------------------------------------------------------ the pair

/// The meeting of two cylinders in turned frames whose axes cross (`x` of
/// operand 0, `y` of operand 1).
pub(super) fn crossing(x: Cyl, y: Cyl, res: f64) -> Result<CylPair> {
    let cyl = [x, y];
    let others = [other_of(y.0, y.1, y.2), other_of(x.0, x.1, x.2)];
    let disc = [
        discriminant(cyl[0], &others[0]),
        discriminant(cyl[1], &others[1]),
    ];
    let piece = |k: usize, plus: bool, range: Option<[[Qd; 2]; 2]>| {
        MeetCrv::new(k, cyl[k], &others[k], plus, range)
    };
    // Charts with a negative antipode, each discriminant's roots, and the
    // near-node test.
    let mut charts: Vec<Option<Chart>> = Vec::new();
    for (a, d) in &disc {
        let chart = negative_chart(d)?;
        let test = chart.clone().unwrap_or(Chart {
            c0: int(1),
            s0: zero(),
        });
        near_node(a, d, &test, res)?;
        charts.push(chart);
    }
    // D_K >= 0 all round (no negative point): rings over K.
    for (k, chart) in charts.iter().enumerate() {
        if chart.is_none() {
            return Ok(CylPair::Quartic(Box::new(Quartic {
                pieces: vec![piece(k, true, None), piece(k, false, None)],
                switches: Vec::new(),
                perpendicular: false,
            })));
        }
    }
    let chart_a = charts[0].clone().expect("a chart");
    let chart_b = charts[1].clone().expect("a chart");
    let pa = disc[0].1.poly(&chart_a);
    let pb = disc[1].1.poly(&chart_b);
    let (mut ra, mut rb) = (roots(&pa)?, roots(&pb)?);
    // Isolators narrowed far below binary64, so their midpoints place the
    // switches.
    for r in ra.iter_mut().chain(rb.iter_mut()) {
        r.refine_for_signs(160);
    }
    if ra.is_empty() || rb.is_empty() {
        // D < 0 all round for one cylinder: no ruling meets the other.
        return Ok(CylPair::Apart);
    }
    if ra.len() % 2 != 0 || rb.len() % 2 != 0 {
        return Err(limit("an odd count of a cylinders' turning points"));
    }
    // Turning points of B's graph: where B's ruling touches A, in A's
    // chart with A's branch (binary64 views; verified below).
    let probe_a = piece(0, true, None);
    let probe_b = piece(1, true, None);
    let mut b_turns: Vec<(f64, bool)> = Vec::new();
    for r in &rb {
        let cs = chart_b.at(&middle(r));
        let (f, c, rr) = cyl[1];
        let base = f.point(&(&c[0] + rr * &cs[0]), &(&c[1] + rr * &cs[1]), &zero());
        // The double root's height: w = -B / A.
        let o = &others[1];
        let n: Vec<R> = (0..2).map(|i| dot(&o.g[i], &f.n)).collect();
        let s: Vec<R> = (0..2).map(|i| dot(&o.g[i], &base) - &o.e[i]).collect();
        let a2 = &n[0] * &n[0] + &n[1] * &n[1];
        let w = -(&n[0] * &s[0] + &n[1] * &s[1]) / a2;
        let pt = qv(&add(&base, &scale(&f.n, &w)));
        let place = probe_a.place(&pt);
        let t = chart_a
            .t_of(&place)
            .ok_or(limit("a turning point at a chart's antipode"))?;
        let branch = match probe_a.branch_sign(&pt) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => return Err(limit("a turning point of both graphs")),
        };
        b_turns.push((t.to_f64(), branch));
    }
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for w in ra.chunks(2) {
        let (a, b) = (&w[0], &w[1]);
        let (af, bf) = (rational_f64(&middle(a)), rational_f64(&middle(b)));
        // Events along the loop: lambda in [0, 2), + branch from a to b,
        // then - branch back.
        let mut events: Vec<(f64, bool)> = vec![(0.0, true), (1.0, true)];
        for &(t, plus) in &b_turns {
            if t <= af || t >= bf {
                continue;
            }
            let l = (t - af) / (bf - af);
            events.push((if plus { l } else { 2.0 - l }, false));
        }
        events.sort_by(|x, y| x.0.total_cmp(&y.0));
        if events.len() < 3 {
            return Err(limit("a loop without turning points of the other graph"));
        }
        let n = events.len();
        for i in 0..n {
            let (l0, l1) = (
                events[i].0,
                events[(i + 1) % n].0 + if i + 1 == n { 2.0 } else { 0.0 },
            );
            if l1 - l0 < 1e-9 {
                return Err(Error::Degenerate("two turning points within rounding"));
            }
        }
        // (lambda) -> (t, branch).
        let at = |l: f64| -> (R, bool) {
            let l = l.rem_euclid(2.0);
            if l < 1.0 {
                (q(af + l * (bf - af)), true)
            } else {
                (q(bf - (l - 1.0) * (bf - af)), false)
            }
        };
        let point = |l: f64| -> Result<(QV, R, bool)> {
            let (t, plus) = at(l);
            if a.compare_rational(&t) != Ordering::Less
                || b.compare_rational(&t) != Ordering::Greater
            {
                return Err(limit("a switch point off its loop"));
            }
            let p = piece(0, plus, None)
                .at(&chart_a.at(&t))
                .ok_or(limit("a switch point off its piece"))?;
            Ok((p, t, plus))
        };
        // Switches between adjacent events of different kinds.
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
            return Err(limit("a loop with one switch"));
        }
        for i in 0..m {
            let (s0, s1) = (&sw[i], &sw[(i + 1) % m]);
            let l1 = s1.0 + if i + 1 == m { 2.0 } else { 0.0 };
            // The run's events.
            let inner: Vec<&(f64, bool)> = events
                .iter()
                .filter(|e| {
                    let l = if e.0 < s0.0 { e.0 + 2.0 } else { e.0 };
                    l > s0.0 && l < l1
                })
                .collect();
            let Some(first) = inner.first() else {
                return Err(limit("a run without turning points"));
            };
            if first.1 {
                // A's turning points: a graph over B's angle.
                let lm = 0.5
                    * (s0.0
                        + if first.0 < s0.0 {
                            first.0 + 2.0
                        } else {
                            first.0
                        });
                let (mid, _, _) = point(lm)?;
                let sign = probe_b.branch_sign(&mid);
                if sign == Ordering::Equal
                    || probe_b.branch_sign(&s0.1) != sign
                    || probe_b.branch_sign(&s1.1) != sign
                {
                    return Err(limit("a run over B's angle changing branch"));
                }
                let (p0, p1, pm) = (
                    probe_b.place(&s0.1),
                    probe_b.place(&s1.1),
                    probe_b.place(&mid),
                );
                let range = if between_ccw(&p0, &pm, &p1) {
                    [p0, p1]
                } else {
                    [p1, p0]
                };
                verify(&pb, &chart_b, &range)?;
                pieces.push(piece(1, sign == Ordering::Greater, Some(range)));
            } else {
                // B's turning points: a graph over A's angle, one branch.
                if s0.3 != s1.3 {
                    return Err(limit("a run over A's angle changing branch"));
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
                    return Err(limit("a turning point inside a piece"));
                }
                let rng = [lo, hi].map(|t| {
                    let cs = chart_a.at(t);
                    [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]
                });
                pieces.push(piece(0, s0.3, Some(rng)));
            }
        }
        switches.extend(sw.into_iter().map(|s| s.1));
    }
    Ok(CylPair::Quartic(Box::new(Quartic {
        pieces,
        switches,
        perpendicular: false,
    })))
}

/// No root of a carrier's discriminant within a counter-clockwise range
/// (the range clear of the chart's antipode): exact Sturm counts at its
/// ends.
fn verify(p: &Poly, chart: &Chart, range: &[[Qd; 2]; 2]) -> Result<()> {
    let (Some(t0), Some(t1)) = (chart.t_of(&range[0]), chart.t_of(&range[1])) else {
        return Err(limit("a piece through a chart's antipode"));
    };
    if t0.cmp(&t1) != Ordering::Less {
        return Err(limit("a piece through a chart's antipode"));
    }
    let chain = sturm(&trim(p.clone()));
    if changes(&chain, &t0) != changes(&chain, &t1) {
        return Err(limit("a turning point inside a piece"));
    }
    Ok(())
}

// ------------------------------------------------------------ cap circles

/// A root of a cap's circle against the other cylinder, in the circle's
/// chart about `(1, 0)`, or its antipode.
#[derive(Debug, Clone)]
pub(super) enum Hit {
    Root(AlgebraicRoot),
    Antipode,
}

/// Where the circle `c + a cos + b sin` meets a cylinder (in turned
/// frames): the roots of its quartic.
pub(super) fn circle_hits(c: &V, a: &V, b: &V, f: &Affine, cy: &P2, ry: &R) -> Result<Vec<Hit>> {
    let o = other_of(f, cy, ry);
    let lin: Vec<Lin> = (0..2)
        .map(|i| [dot(&o.g[i], c) - &o.e[i], dot(&o.g[i], a), dot(&o.g[i], b)])
        .collect();
    let mut form = square_sum(&lin);
    form.k -= &o.r * &o.r;
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let mut p = form.poly(&chart);
    if p.is_empty() {
        return Err(tangency());
    }
    let mut out = Vec::new();
    if form.value(&[int(-1), zero()]) == zero() {
        out.push(Hit::Antipode);
        p = trim(p);
    }
    out.extend(roots(&p)?.into_iter().map(Hit::Root));
    Ok(out)
}

/// Whether a hit lies on an arc from `p` to `q` (directions, rational),
/// counter-clockwise or not, ends included.
pub(super) fn on_arc(hit: &Hit, p: &[R; 2], q2: &[R; 2], ccw: bool) -> bool {
    let (p, q2) = if ccw { (p, q2) } else { (q2, p) };
    let t = |d: &[R; 2]| {
        let den = int(1) + &d[0];
        (den != zero()).then(|| &d[1] / den)
    };
    let (tp, tq) = (t(p), t(q2));
    let cmp = |r: &AlgebraicRoot, x: &R| r.compare_rational(x);
    match hit {
        Hit::Antipode => match (&tp, &tq) {
            (None, _) | (_, None) => true,
            (Some(a), Some(b)) => a > b,
        },
        Hit::Root(r) => match (&tp, &tq) {
            (None, Some(b)) => cmp(r, b) != Ordering::Greater,
            (Some(a), None) => cmp(r, a) != Ordering::Less,
            (None, None) => true,
            (Some(a), Some(b)) => {
                let (ge, le) = (cmp(r, a) != Ordering::Less, cmp(r, b) != Ordering::Greater);
                if a <= b {
                    ge && le
                } else {
                    ge || le
                }
            }
        },
    }
}
