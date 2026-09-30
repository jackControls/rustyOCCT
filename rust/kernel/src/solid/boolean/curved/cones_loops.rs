//! S9d.3c: two ruled faces meeting in loops (REVIEW_NOTES.md, S9d.3c
//! refined): a cone against a cylinder or another cone whose rulings'
//! discriminant `D = B^2 - A C` changes sign over both carriers' angles, in
//! any affine frames. As S9c.2b.1 meets two turned cylinders
//! (`turned::crossing`), generalised to ruled carriers: each interval of
//! `D_0 > 0` over carrier 0's angle is one component (its plus branch up,
//! its minus branch back), the other carrier's turning points placed along
//! it by binary64 views, switches at rational angles of carrier 0 between
//! adjacent turning points of different kinds; the pieces between them are
//! graphs over one carrier's angle (`Curve3::Meet`), each verified exactly
//! (no root of its carrier's `D` or `A` inside its range, its ends and an
//! exact interior point on its branch, the boundaries ordered exactly along
//! the component). Where the rulings reach the other's asymptotic
//! directions (`A`'s simple real roots: two cones whose direction cones
//! cross) a component runs through infinity: it is cut there on the branch
//! that runs off, the finite branch's point a switch, as S9d.3b.2's open
//! pieces.
use super::graph::between_ccw;
use super::meet::{tangency, CylPair};
use super::num::*;
use super::procedural::{MeetCrv, Other, Quartic, Ruled};
use super::turned::{
    changes, middle, near_node, negative_chart, roots, ruled_quadratic, sturm, trim, Chart, Form,
    Poly,
};
use crate::polynomial::real::AlgebraicRoot;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

/// A root of a carrier's `A` (a direction along the other's asymptotic
/// cone): its chart `t` and direction (in `Q(alpha)`), whether the plus
/// branch runs off to infinity there (`-sign B`), and the other branch's
/// point (`w = -C / 2B`).
struct Asym {
    t: Qd,
    t64: f64,
    cs: [Qd; 2],
    inf_plus: bool,
    finite: QV,
}

/// One of the meeting's two carriers.
struct Carrier<'a> {
    op: usize,
    ruled: Ruled<'a>,
    other: Other,
    a: Form,
    b: Form,
    chart: Chart,
    /// `D` and `A` in the chart (`A` empty when constant), `D`'s real roots.
    dp: Poly,
    ap: Poly,
    droots: Vec<AlgebraicRoot>,
    asym: Vec<Asym>,
}

/// `A`'s size over the turn, sampled at rational directions: its least
/// (halved) where it has no real root, else its largest; the near-node
/// test's scale (a guard, not a decision).
fn a_scale(a: &Form, rooted: bool) -> R {
    if let Some(k) = a.as_constant() {
        return if k < zero() { -k } else { k };
    }
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let vals: Vec<f64> = (-16..=16)
        .map(|i| rational_f64(&a.value(&chart.at(&(int(i) / int(4))))).abs())
        .collect();
    if rooted {
        q(vals.iter().fold(0.0f64, |m, &v| m.max(v)))
    } else {
        q(vals.iter().fold(f64::INFINITY, |m, &v| m.min(v)) / 2.0)
    }
}

/// A chart's direction at `t` of any field.
fn chart_at(chart: &Chart, t: &K) -> Result<[Qd; 2]> {
    let den = t.mul(t).add(&K::Rat(int(1)));
    let inv = den.recip().ok_or(limit("a direction at infinity"))?;
    let c = K::Rat(int(1)).sub(&t.mul(t)).mul(&inv);
    let s = t.scale(&int(2)).mul(&inv);
    Ok([
        Qd::of(c.scale(&chart.c0).sub(&s.scale(&chart.s0))),
        Qd::of(c.scale(&chart.s0).add(&s.scale(&chart.c0))),
    ])
}

/// A ruling's point at `w` over a direction of any field.
pub(super) fn ruling_point(ruled: Ruled, cs: &[Qd; 2], w: &Qd) -> QV {
    let f = ruled.f;
    let base = qadd(
        &qv(&f.point(&ruled.c[0], &ruled.c[1], &zero())),
        &qadd(
            &qscale(&f.x, &cs[0].scale(ruled.r)),
            &qscale(&f.y, &cs[1].scale(ruled.r)),
        ),
    );
    let dir = qadd(
        &qv(&f.n),
        &qadd(
            &qscale(&f.x, &cs[0].scale(ruled.k)),
            &qscale(&f.y, &cs[1].scale(ruled.k)),
        ),
    );
    qadd(&base, &dir.map(|x| x.mul(w)))
}

/// A ruling's point at a rational height over a rational direction.
fn ruling_at(ruled: Ruled, cs: &[R; 2], w: &R) -> V {
    let f = ruled.f;
    let radial = add(&scale(&f.x, &cs[0]), &scale(&f.y, &cs[1]));
    let base = add(
        &f.point(&ruled.c[0], &ruled.c[1], &zero()),
        &scale(&radial, ruled.r),
    );
    let dir = add(&f.n, &scale(&radial, ruled.k));
    add(&base, &scale(&dir, w))
}

/// A ruling's direction in binary64 over a direction of any field.
fn ruling_dir64(ruled: Ruled, cs: &[Qd; 2]) -> [f64; 3] {
    let f = ruled.f;
    let (c, s, k) = (cs[0].to_f64(), cs[1].to_f64(), rational_f64(ruled.k));
    [0, 1, 2].map(|j| {
        rational_f64(&f.n[j]) + k * (c * rational_f64(&f.x[j]) + s * rational_f64(&f.y[j]))
    })
}

/// A polynomial's exact value at a number of any field.
fn eval(p: &Poly, x: &Qd) -> Qd {
    let mut acc = Qd::rat(zero());
    for c in p.iter().rev() {
        acc = acc.mul(x).add_r(c);
    }
    acc
}

/// The count of a polynomial's roots in `(t0, t1)`, `t0 < t1` (Sturm's, a
/// root at `t1` taken out).
fn inside(p: &Poly, t0: &Qd, t1: &Qd) -> usize {
    let chain = sturm(&trim(p.clone()));
    let n = changes(&chain, t0).saturating_sub(changes(&chain, t1));
    if eval(p, t1).sign() == Ordering::Equal {
        n.saturating_sub(1)
    } else {
        n
    }
}

impl<'a> Carrier<'a> {
    fn new(op: usize, ruled: Ruled<'a>, other: Other, res: f64) -> Result<Self> {
        let (a, b, c) = ruled_quadratic(ruled, &other);
        let d = b.mul(&b).sub(&a.mul(&c));
        let chart = negative_chart(&d)?.ok_or(limit(
            "a loop's carrier with its discriminant nowhere negative",
        ))?;
        let dp = trim(d.poly(&chart));
        let varies = a.as_constant().is_none();
        let ap = if varies {
            trim(a.poly(&chart))
        } else {
            Vec::new()
        };
        let mut asym = Vec::new();
        if varies {
            let found = match roots(&ap) {
                Ok(r) => r,
                // A repeated direction (two direction cones tangent).
                Err(Error::Degenerate(_)) => {
                    return Err(Error::Degenerate(
                        "a ruling along the other's asymptotic direction",
                    ))
                }
                Err(e) => return Err(e),
            };
            for root in found {
                let g = Arc::new(Gen::new(ap.clone(), root));
                let t = K::generator(&g);
                let t64 = Qd::of(t.clone()).to_f64();
                let cs = chart_at(&chart, &t)?;
                let bv = b.value_q(&cs);
                let sb = bv.sign();
                if sb == Ordering::Equal {
                    // `D = B^2 = 0`: a ruling missing the other along its
                    // asymptotic direction, a tangency at infinity.
                    return Err(tangency());
                }
                let w = bv
                    .scale(&int(2))
                    .recip()
                    .ok_or(tangency())?
                    .mul(&c.value_q(&cs).neg());
                let finite = ruling_point(ruled, &cs, &w);
                asym.push(Asym {
                    t: Qd::of(t),
                    t64,
                    cs,
                    inf_plus: sb == Ordering::Less,
                    finite,
                });
            }
        }
        near_node(&a_scale(&a, !asym.is_empty()), &d, &chart, res)?;
        let mut droots = roots(&dp)?;
        for r in droots.iter_mut() {
            r.refine_for_signs(160);
        }
        Ok(Self {
            op,
            ruled,
            other,
            a,
            b,
            chart,
            dp,
            ap,
            droots,
            asym,
        })
    }

    fn piece(&self, plus: bool, range: Option<[[Qd; 2]; 2]>) -> MeetCrv {
        MeetCrv::ruled(self.op, self.ruled, &self.other, plus, range)
    }

    /// No root of `D` or `A` strictly inside a counter-clockwise range (clear
    /// of the chart's antipode, where `D < 0`).
    fn verify(&self, range: &[[Qd; 2]; 2]) -> Result<()> {
        let (Some(t0), Some(t1)) = (self.chart.t_of(&range[0]), self.chart.t_of(&range[1])) else {
            return Err(limit("a piece through a chart's antipode"));
        };
        if t0.cmp(&t1) != Ordering::Less {
            return Err(limit("a piece through a chart's antipode"));
        }
        if inside(&self.dp, &t0, &t1) != 0 {
            return Err(limit("a turning point inside a piece"));
        }
        if !self.ap.is_empty() && inside(&self.ap, &t0, &t1) != 0 {
            return Err(limit("an asymptotic direction inside a piece"));
        }
        Ok(())
    }
}

/// What bounds a run along a component.
#[derive(Clone)]
enum End {
    /// A switch at a rational angle of carrier 0.
    Switch(Box<QV>),
    /// Carrier 0's asymptotic direction `i`: its branch running off.
    Cut(usize),
    /// Carrier 0's asymptotic direction `i`: its finite branch's point.
    Fin0(usize),
    /// Carrier 1's asymptotic direction `i`: its finite branch's point.
    Fin1(usize),
}

/// A boundary: its binary64 place `l` along the component, its exact key
/// (branch, chart `t` on carrier 0) and what it is.
#[derive(Clone)]
struct Bd {
    l: f64,
    key: (bool, Qd),
    end: End,
}

/// The exact order along a component: the plus branch up, then the minus
/// branch back.
fn key_cmp(x: &(bool, Qd), y: &(bool, Qd)) -> Ordering {
    match (x.0, y.0) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        (true, true) => x.1.cmp(&y.1),
        (false, false) => y.1.cmp(&x.1),
    }
}

/// A turning point along a component: of carrier 0 (`true`) or 1.
#[derive(Clone, Copy)]
struct Turn {
    l: f64,
    zero: bool,
}

/// Two ruled carriers (`first` carrier 0, a cylinder where there is one)
/// meeting in loops, or in components through infinity: their pieces and
/// switches.
pub(super) fn loops(
    first: (usize, Ruled, Other),
    second: (usize, Ruled, Other),
    res: f64,
) -> Result<CylPair> {
    let car = [
        Carrier::new(first.0, first.1, first.2, res)?,
        Carrier::new(second.0, second.1, second.2, res)?,
    ];
    if car[0].droots.is_empty() || car[1].droots.is_empty() {
        // D < 0 all round for one carrier: no ruling meets the other.
        return Ok(CylPair::Apart);
    }
    if car[0].droots.len() % 2 != 0 || car[1].droots.len() % 2 != 0 {
        return Err(limit("an odd count of a carrier's turning points"));
    }
    let probe = [car[0].piece(true, None), car[1].piece(true, None)];
    let branch = |k: usize, p: &QV| -> Result<bool> {
        match probe[k].branch_sign(p) {
            Ordering::Greater => Ok(true),
            Ordering::Less => Ok(false),
            Ordering::Equal => Err(limit("a turning point of both graphs")),
        }
    };
    // Carrier 1's turning points on carrier 0: its ruling touching the first
    // quadric (`w = -B / A`), binary64 views of chart `t` and branch.
    let mut turns1: Vec<(f64, bool)> = Vec::new();
    for r in &car[1].droots {
        let cs = car[1].chart.at(&middle(r));
        let (av, bv) = (car[1].a.value(&cs), car[1].b.value(&cs));
        if av == zero() {
            return Err(limit("a turning point along an asymptotic direction"));
        }
        let pt = qv(&ruling_at(car[1].ruled, &cs, &(-bv / av)));
        let t = car[0]
            .chart
            .t_of(&probe[0].place(&pt))
            .ok_or(limit("a turning point at a chart's antipode"))?;
        turns1.push((t.to_f64(), branch(0, &pt)?));
    }
    // Carrier 1's finite asymptotic points on carrier 0, exactly.
    let mut fin1: Vec<(Qd, bool)> = Vec::new();
    for s in &car[1].asym {
        let t = car[0]
            .chart
            .t_of(&probe[0].place(&s.finite))
            .ok_or(limit("an asymptotic point at a chart's antipode"))?;
        fin1.push((t, branch(0, &s.finite)?));
    }
    let chain0 = sturm(&car[0].dp);
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for w in car[0].droots.chunks(2) {
        let (a, b) = (&w[0], &w[1]);
        let (af, bf) = (rational_f64(&middle(a)), rational_f64(&middle(b)));
        let within = |t: f64| t > af && t < bf;
        let lam = |t: f64, plus: bool| {
            let l = (t - af) / (bf - af);
            if plus {
                l
            } else {
                2.0 - l
            }
        };
        // A rational point inside the interval, and whether an exact chart
        // `t` lies in it (no root of `D_0` between them).
        let tref = (a.isolator().1 + b.isolator().0) / int(2);
        let tref = Qd::rat(tref);
        let same = |t: &Qd| -> bool {
            let (lo, hi) = if t.cmp(&tref) == Ordering::Less {
                (t, &tref)
            } else {
                (&tref, t)
            };
            changes(&chain0, lo) == changes(&chain0, hi)
        };
        let mut turns: Vec<Turn> = vec![Turn { l: 0.0, zero: true }, Turn { l: 1.0, zero: true }];
        for &(t, plus) in &turns1 {
            if within(t) {
                turns.push(Turn {
                    l: lam(t, plus),
                    zero: false,
                });
            }
        }
        let mut bd: Vec<Bd> = Vec::new();
        for (i, s) in car[0].asym.iter().enumerate() {
            if !within(s.t64) {
                continue;
            }
            bd.push(Bd {
                l: lam(s.t64, s.inf_plus),
                key: (s.inf_plus, s.t.clone()),
                end: End::Cut(i),
            });
            bd.push(Bd {
                l: lam(s.t64, !s.inf_plus),
                key: (!s.inf_plus, s.t.clone()),
                end: End::Fin0(i),
            });
        }
        for (i, (t, plus)) in fin1.iter().enumerate() {
            let tf = t.to_f64();
            if within(tf) {
                bd.push(Bd {
                    l: lam(tf, *plus),
                    key: (*plus, t.clone()),
                    end: End::Fin1(i),
                });
            }
        }
        // Every event apart from the others by more than rounding.
        let mut ls: Vec<f64> = turns
            .iter()
            .map(|x| x.l)
            .chain(bd.iter().map(|x| x.l))
            .collect();
        ls.sort_by(f64::total_cmp);
        let n = ls.len();
        for i in 0..n {
            let l1 = ls[(i + 1) % n] + if i + 1 == n { 2.0 } else { 0.0 };
            if l1 - ls[i] < 1e-9 {
                return Err(Error::Degenerate("two turning points within rounding"));
            }
        }
        // (lambda) -> (t, branch), a rational point of carrier 0's chart.
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
            let p = car[0]
                .piece(plus, None)
                .at(&car[0].chart.at(&t))
                .ok_or(limit("a switch point off its piece"))?;
            Ok((p, t, plus))
        };
        // Switches between adjacent turning points of different kinds with
        // no boundary between them.
        enum Item {
            T(Turn),
            B,
        }
        let mut items: Vec<(f64, Item)> = turns
            .iter()
            .map(|t| (t.l, Item::T(*t)))
            .chain(bd.iter().map(|x| (x.l, Item::B)))
            .collect();
        items.sort_by(|x, y| x.0.total_cmp(&y.0));
        let n = items.len();
        for i in 0..n {
            let (Item::T(t0), Item::T(t1)) = (&items[i].1, &items[(i + 1) % n].1) else {
                continue;
            };
            if t0.zero == t1.zero {
                continue;
            }
            let l1 = t1.l + if i + 1 == n { 2.0 } else { 0.0 };
            let l = 0.5 * (t0.l + l1);
            let (p, t, plus) = point(l)?;
            bd.push(Bd {
                l: l.rem_euclid(2.0),
                key: (plus, Qd::rat(t)),
                end: End::Switch(Box::new(p)),
            });
        }
        if turns.iter().all(|t| t.zero) && bd.is_empty() {
            return Err(limit("a loop without turning points of the other graph"));
        }
        // The boundaries in their exact order along the component, each of
        // it (no root of `D_0` between it and the interval's point).
        bd.sort_by(|x, y| key_cmp(&x.key, &y.key));
        if bd.iter().any(|x| !same(&x.key.1)) {
            return Err(limit("a boundary off its loop"));
        }
        let m = bd.len();
        if m < 2 {
            return Err(limit("a component with one boundary"));
        }
        for w in bd.windows(2) {
            if key_cmp(&w[0].key, &w[1].key) != Ordering::Less {
                return Err(limit("two boundaries at one point"));
            }
        }
        for i in 0..m {
            let (s0, s1) = (&bd[i], &bd[(i + 1) % m]);
            let wrap = i + 1 == m;
            let l1 = s1.l + if wrap { 2.0 } else { 0.0 };
            let run: Vec<&Turn> = turns
                .iter()
                .filter(|t| {
                    let l = if t.l < s0.l { t.l + 2.0 } else { t.l };
                    l > s0.l && l < l1
                })
                .collect();
            let zero_kind = run.iter().any(|t| t.zero);
            if zero_kind && run.iter().any(|t| !t.zero) {
                return Err(limit("a run with turning points of both kinds"));
            }
            if zero_kind {
                // Carrier 0's turning points: a graph over carrier 1's angle,
                // through an exact point of the run.
                let first = run
                    .iter()
                    .map(|t| if t.l < s0.l { t.l + 2.0 } else { t.l })
                    .fold(f64::INFINITY, f64::min);
                let (mid, t, plus) = point(0.5 * (s0.l + first))?;
                let key = (plus, Qd::rat(t));
                let between = if wrap {
                    key_cmp(&s0.key, &key) == Ordering::Less
                        || key_cmp(&key, &s1.key) == Ordering::Less
                } else {
                    key_cmp(&s0.key, &key) == Ordering::Less
                        && key_cmp(&key, &s1.key) == Ordering::Less
                };
                if !between {
                    return Err(limit("a run's point off the run"));
                }
                let plus1 = branch(1, &mid)?;
                let e0 = end_on_one(&car, &probe, s0, plus1)?;
                let e1 = end_on_one(&car, &probe, s1, plus1)?;
                let pm = probe[1].place(&mid);
                let range = if between_ccw(&e0, &pm, &e1) {
                    [e0, e1]
                } else {
                    [e1, e0]
                };
                car[1].verify(&range)?;
                pieces.push(car[1].piece(plus1, Some(range)));
            } else {
                // Carrier 1's turning points (or none): a graph over carrier
                // 0's angle on one branch.
                if s0.key.0 != s1.key.0 {
                    return Err(limit("a run over carrier 0's angle changing branch"));
                }
                let plus = s0.key.0;
                let (d0, d1) = (end_on_zero(&car, &probe, s0), end_on_zero(&car, &probe, s1));
                let range = if plus { [d0, d1] } else { [d1, d0] };
                car[0].verify(&range)?;
                pieces.push(car[0].piece(plus, Some(range)));
            }
        }
        for s in bd {
            match s.end {
                End::Switch(p) => switches.push(*p),
                End::Fin0(i) => switches.push(car[0].asym[i].finite.clone()),
                End::Fin1(i) => switches.push(car[1].asym[i].finite.clone()),
                End::Cut(_) => {}
            }
        }
    }
    Ok(CylPair::Quartic(Box::new(Quartic { pieces, switches })))
}

/// A boundary's direction on carrier 0.
fn end_on_zero(car: &[Carrier; 2], probe: &[MeetCrv; 2], s: &Bd) -> [Qd; 2] {
    match &s.end {
        End::Switch(p) => probe[0].place(p),
        End::Cut(i) | End::Fin0(i) => car[0].asym[*i].cs.clone(),
        End::Fin1(i) => probe[0].place(&car[1].asym[*i].finite),
    }
}

/// A boundary's direction on carrier 1, on the branch `plus`: a point's
/// place there (on that branch), carrier 1's own asymptotic direction (its
/// finite branch that one), or for a cut the asymptotic direction of carrier
/// 1 along the same line at infinity (its branch running off that one; the
/// verified piece reaching infinity there only, the first point at infinity
/// along the run, whichever direction is taken).
fn end_on_one(car: &[Carrier; 2], probe: &[MeetCrv; 2], s: &Bd, plus: bool) -> Result<[Qd; 2]> {
    let on = |p: &QV| -> Result<[Qd; 2]> {
        let want = if plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if probe[1].branch_sign(p) != want {
            return Err(limit("a run over carrier 1's angle changing branch"));
        }
        Ok(probe[1].place(p))
    };
    match &s.end {
        End::Switch(p) => on(p),
        End::Fin0(i) => on(&car[0].asym[*i].finite),
        End::Fin1(i) => {
            let x = &car[1].asym[*i];
            if x.inf_plus == plus {
                return Err(limit("a run over carrier 1's angle changing branch"));
            }
            Ok(x.cs.clone())
        }
        End::Cut(i) => {
            let d0 = ruling_dir64(car[0].ruled, &car[0].asym[*i].cs);
            let unit = |d: [f64; 3]| {
                let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                d.map(|x| x / l)
            };
            let d0 = unit(d0);
            let off = |d: [f64; 3]| {
                let d = unit(d);
                let c = [
                    d0[1] * d[2] - d0[2] * d[1],
                    d0[2] * d[0] - d0[0] * d[2],
                    d0[0] * d[1] - d0[1] * d[0],
                ];
                c[0] * c[0] + c[1] * c[1] + c[2] * c[2]
            };
            let x = car[1]
                .asym
                .iter()
                .min_by(|x, y| {
                    off(ruling_dir64(car[1].ruled, &x.cs))
                        .total_cmp(&off(ruling_dir64(car[1].ruled, &y.cs)))
                })
                .ok_or(limit("a cut without the other's asymptotic direction"))?;
            if x.inf_plus != plus {
                return Err(limit("a run over carrier 1's angle changing branch"));
            }
            Ok(x.cs.clone())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sturm's count of the roots strictly inside a range, a root at either
    /// end left out (a piece ending at an asymptotic direction).
    #[test]
    fn roots_inside_a_range_leave_its_ends_out() {
        // (t - 1)(t - 2)(t - 3) = t^3 - 6 t^2 + 11 t - 6.
        let p: Poly = vec![int(-6), int(11), int(-6), int(1)];
        let r = |x: i64, y: i64| Qd::rat(R::new(x.into(), y.into()));
        assert_eq!(inside(&p, &r(0, 1), &r(4, 1)), 3);
        assert_eq!(inside(&p, &r(1, 1), &r(3, 1)), 1);
        assert_eq!(inside(&p, &r(3, 2), &r(2, 1)), 0);
        assert_eq!(inside(&p, &r(1, 1), &r(5, 2)), 1);
        // An end in `Q(sqrt 2)` between the roots.
        let s = Qd::new(zero(), int(1), int(2));
        assert_eq!(inside(&p, &s, &r(3, 1)), 1);
    }

    /// Two boundaries' order along a component: the plus branch up, then
    /// the minus branch back.
    #[test]
    fn boundaries_order_along_a_component() {
        let k = |plus: bool, t: i64| (plus, Qd::rat(int(t)));
        assert_eq!(key_cmp(&k(true, 1), &k(true, 2)), Ordering::Less);
        assert_eq!(key_cmp(&k(false, 1), &k(false, 2)), Ordering::Greater);
        assert_eq!(key_cmp(&k(true, 5), &k(false, -5)), Ordering::Less);
    }
}
