//! Certified evaluation of a spline wall's meeting with a cylinder
//! (`Curve3::WallMeet`, S9f.2b; D13), in either tier and over any of the
//! integrands' numbers (`quadrature::Num`: plain enclosures, Taylor series,
//! jets).
//!
//! The wall is its face's stored B-spline surface of degree one in `v`: on
//! each knot span its two pole rows are exact Bézier rows (the surface's
//! exact patches), so its ruling there is `L(ū) + t M(ū)` with `L` the low
//! row's Bernstein polynomial, `M` the rows' difference over the `v` range
//! and `ū` the span's local parameter. The other cylinder's function along
//! the ruling is `a t^2 + 2 b t + c`, and the curve's `t` is `(-b + s
//! sqrt(b^2 - a c)) / a` or `c / (-b - s sqrt(b^2 - a c))`, the one whose
//! middle cancels less (the binary64 curve's own choice, both the same
//! value). `a`, `b`, `c` and the discriminant `d = a r^2 - (P x M)^2` (`P`
//! and `M` the ruling's foot and direction in the cylinder's axes:
//! Lagrange's identity) are polynomials in `ū` whose Bernstein coefficients
//! are made exactly once per wall and cylinder (products of Bernstein
//! polynomials), so each is enclosed at a point within a few units in the
//! last place of its coefficients: near a turning point, where `v` is `d`'s
//! square root, `d` formed from enclosed factors (`P` and `M` each from two
//! rows) was a hundred times wider, and its root's enclosure as much. Over
//! a range (a remainder's box) they are formed from the factors' exact
//! polynomials instead, which overestimate less than the products' higher
//! degrees. A knot span's polynomial holds the curve
//! on that span only: across a knot the wall is C^k (k its least interior
//! continuity, at least one by R4), so a jet over a base across a knot is
//! the union of both spans' jets up to order `k + 1` (a Taylor remainder's
//! bound holds for a function whose `k`-th derivative is absolutely
//! continuous) and `None` past it; integrals along the curve are split at
//! the knots exactly (`pieces`), each piece on its span's polynomial.
use super::quadrature::Num;
use crate::certified::{Fast, Real};
use crate::jet::Jet;
use crate::topology::WallMeet;
use crate::BSplineSurface3;
use num_rational::BigRational as R;

use std::cell::RefCell;
use std::sync::Arc;

/// Exact Bernstein coefficients and their binary64 enclosures, made once
/// (the binary64 tier's conversions are most of an evaluation otherwise).
struct Coefficients {
    exact: Vec<R>,
    fast: Vec<Fast>,
}

impl Coefficients {
    fn new(exact: Vec<R>) -> Self {
        let fast = exact.iter().map(Fast::from_r).collect();
        Self { exact, fast }
    }

    fn at<T: Real, N: Num<T>>(&self, t: &N) -> N {
        let b: Vec<T> = self
            .exact
            .iter()
            .zip(&self.fast)
            .map(|(x, f)| cached::<T>(f, x))
            .collect();
        bernstein(&b, t)
    }
}

/// A knot span of the wall met by one cylinder: its `u` range (exact and
/// binary64), the wall's `v` range, the low row's and the rulings'
/// direction's Bernstein coordinates (`low[k]`, `dir[k]`: coordinate `k`)
/// and the polynomials `a`, `b`, `c` and `d` in `ū`.
pub(super) struct Span {
    pub(super) u: [R; 2],
    pub(super) uf: [f64; 2],
    pub(super) v: [R; 2],
    low: [Coefficients; 3],
    dir: [Coefficients; 3],
    abcd: [Coefficients; 4],
    /// The foot's and the direction's coordinates along the cylinder's
    /// `x` and `y` (`wx`, `wy`, `mx`, `my`), for ranges.
    axes: [Coefficients; 4],
    /// `-u0`, `1 / (u1 - u0)` and `v0`, exact, and binary64.
    constants: [R; 3],
    fast_constants: [Fast; 3],
}

/// The wall's spans met by the cylinder and the wall's least continuity
/// across an interior knot.
pub(super) struct Spans {
    pub(super) spans: Vec<Span>,
    continuity: usize,
}

/// Walls and cylinders kept (each pair's polynomials are exact rational
/// work).
const KEPT: usize = 64;

type Key = (BSplineSurface3, crate::Frame3, u64);

thread_local! {
    static WALLS: RefCell<Vec<(Key, Arc<Spans>)>> = const { RefCell::new(Vec::new()) };
}

/// The coefficients of the product of two Bernstein polynomials.
fn product(f: &[R], g: &[R]) -> Vec<R> {
    let (m, n) = (f.len() - 1, g.len() - 1);
    let binomial = |n: usize, k: usize| -> R {
        let mut out = one();
        for i in 0..k {
            out = out * R::from_integer((n - i).into()) / R::from_integer((i + 1).into());
        }
        out
    };
    (0..=m + n)
        .map(|k| {
            let mut sum = zero();
            for i in k.saturating_sub(n)..=k.min(m) {
                sum += binomial(m, i) * binomial(n, k - i) * &f[i] * &g[k - i];
            }
            sum / binomial(m + n, k)
        })
        .collect()
}

fn combine(f: &[R], g: &[R], scale: &R) -> Vec<R> {
    f.iter().zip(g).map(|(a, b)| a + b * scale).collect()
}

/// The wall's spans met by the curve's cylinder, made once per wall and
/// cylinder (`None` for a wall that is not a nonrational surface of degree
/// one in `v` over one `v` span).
pub(super) fn spans(m: &WallMeet) -> Option<Arc<Spans>> {
    let key: Key = (m.wall.clone(), m.other, m.other_radius.to_bits());
    if let Some(hit) = WALLS.with(|w| {
        w.borrow()
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, x)| x.clone())
    }) {
        return Some(hit);
    }
    let exact = |x: f64| R::from_float(x);
    let (o, x2, y2) = (
        m.other.origin().to_array(),
        m.other.x().to_array(),
        m.other.y().to_array(),
    );
    let o: [R; 3] = [exact(o[0])?, exact(o[1])?, exact(o[2])?];
    let x2: [R; 3] = [exact(x2[0])?, exact(x2[1])?, exact(x2[2])?];
    let y2: [R; 3] = [exact(y2[0])?, exact(y2[1])?, exact(y2[2])?];
    let r = exact(m.other_radius)?;
    let r2 = &r * &r;
    let s = &m.wall;
    if s.is_rational() || s.u_knots().is_periodic() || s.v_knots().degree() != 1 {
        return None;
    }
    if s.v_knots().pole_count() != 2 {
        return None;
    }
    let patches = s.bezier_patches().ok()?;
    let mut spans = Vec::new();
    for q in &patches {
        let [du, dv] = q.degrees();
        if dv != 1 {
            return None;
        }
        let poles = q.homogeneous_poles();
        let row = |j: usize| -> [Vec<R>; 3] {
            std::array::from_fn(|k| {
                (0..=du)
                    .map(|i| {
                        let p = &poles[i * 2 + j];
                        &p[k] / &p[3]
                    })
                    .collect()
            })
        };
        let [u, v] = q.domain().clone();
        let f = |x: &R| crate::solid::split::rational_f64(x);
        let (low, high) = (row(0), row(1));
        let per = one() / (&v[1] - &v[0]);
        let dir: [Vec<R>; 3] = std::array::from_fn(|k| {
            low[k]
                .iter()
                .zip(&high[k])
                .map(|(a, b)| (b - a) * &per)
                .collect()
        });
        // The foot's and the direction's coordinates in the cylinder's
        // axes, then `a`, `b`, `c` and `d`.
        let axis = |p: &[Vec<R>; 3], e: &[R; 3], shift: bool| -> Vec<R> {
            (0..=du)
                .map(|i| {
                    (0..3)
                        .map(|k| {
                            let x = if shift {
                                &p[k][i] - &o[k]
                            } else {
                                p[k][i].clone()
                            };
                            x * &e[k]
                        })
                        .fold(zero(), |s, x| s + x)
                })
                .collect()
        };
        let (wx, wy) = (axis(&low, &x2, true), axis(&low, &y2, true));
        let (mx, my) = (axis(&dir, &x2, false), axis(&dir, &y2, false));
        let minus = -one();
        let a = combine(&product(&mx, &mx), &product(&my, &my), &one());
        let b = combine(&product(&wx, &mx), &product(&wy, &my), &one());
        let c: Vec<R> = combine(&product(&wx, &wx), &product(&wy, &wy), &one())
            .into_iter()
            .map(|x| x - &r2)
            .collect();
        let cross = combine(&product(&wx, &my), &product(&wy, &mx), &minus);
        // `a r^2` raised to `cross^2`'s degree (times one of degree
        // `2 du`).
        let ar2: Vec<R> = a.iter().map(|x| x * &r2).collect();
        let d = combine(
            &product(&ar2, &vec![one(); 2 * du + 1]),
            &product(&cross, &cross),
            &minus,
        );
        let constants = [-&u[0], one() / (&u[1] - &u[0]), v[0].clone()];
        let fast_constants = std::array::from_fn(|i| Fast::from_r(&constants[i]));
        spans.push(Span {
            uf: [f(&u[0]), f(&u[1])],
            u,
            v,
            low: low.map(Coefficients::new),
            dir: dir.map(Coefficients::new),
            abcd: [a, b, c, d].map(Coefficients::new),
            axes: [wx, wy, mx, my].map(Coefficients::new),
            constants,
            fast_constants,
        });
    }
    spans.sort_by(|a, b| a.u[0].cmp(&b.u[0]));
    let k = s.u_knots();
    let p = k.degree();
    let m_inner = k.multiplicities()[1..k.multiplicities().len() - 1]
        .iter()
        .copied()
        .max()
        .unwrap_or(0);
    let continuity = p.saturating_sub(m_inner.max(1));
    let out = Arc::new(Spans { spans, continuity });
    WALLS.with(|w| {
        let mut w = w.borrow_mut();
        if w.len() >= KEPT {
            w.drain(..KEPT / 2);
        }
        w.push((key, out.clone()));
    });
    Some(out)
}

/// A Bernstein polynomial at `t` by de Casteljau: at a point as `(1 - t) a +
/// t b` (enclosures that do not grow with the levels where `t` lies in `[0,
/// 1]`; `a + t (b - a)` widens them by up to `1 + 2 t` a level), over a
/// range as `a + t (b - a)` (the other form's two products each take the
/// range's width at the coefficients' size).
pub(super) fn bernstein<T: Real, N: Num<T>>(b: &[T], t: &N) -> N {
    let point = t.sharp();
    let s = t.lift(&T::exact_f64(1.0)).sub(t);
    let mut row: Vec<N> = b.iter().map(|x| t.lift(x)).collect();
    for last in (1..row.len()).rev() {
        for i in 0..last {
            row[i] = if point {
                row[i].mul(&s).add(&row[i + 1].mul(t))
            } else {
                row[i].add(&t.mul(&row[i + 1].sub(&row[i])))
            };
        }
    }
    row.swap_remove(0)
}

/// An exact number in tier `T`: its binary64 enclosure made once where `T`
/// is that tier.
fn cached<T: Real>(fast: &Fast, x: &R) -> T {
    T::of_fast(*fast).unwrap_or_else(|| T::from_r(x))
}

fn fc<T: Real>(x: f64) -> T {
    T::exact_f64(x)
}

/// The wall's `v` and the curve's point at the wall's `u` (any number),
/// on span `span`'s polynomial.
pub(super) fn eval<T: Real, N: Num<T>>(m: &WallMeet, span: &Span, u: &N) -> Option<(N, [N; 3])> {
    let constant = |i: usize| cached::<T>(&span.fast_constants[i], &span.constants[i]);
    let ub = u.shift(&constant(0)).scale(&constant(1));
    let low: [N; 3] = std::array::from_fn(|k| span.low[k].at(&ub));
    let dir: [N; 3] = std::array::from_fn(|k| span.dir[k].at(&ub));
    let [a, b, c, d] = if ub.sharp() {
        std::array::from_fn(|i| span.abcd[i].at::<T, N>(&ub))
    } else {
        // Over a range the polynomials of degree `4p` overestimate (their
        // coefficients far above a small discriminant): the factors'.
        let [wx, wy, mx, my] = std::array::from_fn(|i| span.axes[i].at::<T, N>(&ub));
        let r = fc::<T>(m.other_radius);
        let r2 = r.mul(&r);
        let a = mx.square().add(&my.square());
        let b = wx.mul(&mx).add(&wy.mul(&my));
        let c = wx.square().add(&wy.square()).shift(&r2.neg());
        let cross = wx.mul(&my).sub(&wy.mul(&mx));
        let d = a.scale(&r2).sub(&cross.square());
        [a, b, c, d]
    };
    let sq = d.sqrt()?.scale(&fc(m.sign));
    let (p, q) = (sq.sub(&b), sq.neg().sub(&b));
    let t = if p.mid().abs() >= q.mid().abs() {
        p.div(&a)?
    } else {
        c.div(&q)?
    };
    let point = std::array::from_fn(|k| low[k].add(&dir[k].mul(&t)));
    Some((t.shift(&constant(2)), point))
}

/// The discriminant `d` at a local `ū` in binary64 (no enclosure: a guide
/// to where a range's series can run).
pub(super) fn discriminant_f64(span: &Span, ub: f64) -> f64 {
    let mut row: Vec<f64> = span.abcd[3]
        .fast
        .iter()
        .map(|x| {
            let (lo, hi) = x.bounds_f64();
            0.5 * lo + 0.5 * hi
        })
        .collect();
    for last in (1..row.len()).rev() {
        for i in 0..last {
            row[i] = (1.0 - ub) * row[i] + ub * row[i + 1];
        }
    }
    row[0]
}

/// The spans a `u` enclosure meets (closed).
fn meeting(spans: &Spans, lo: f64, hi: f64) -> Vec<usize> {
    let n = spans.spans.len();
    (0..n)
        .filter(|&k| {
            let [a, b] = spans.spans[k].uf;
            (hi > a && lo < b)
                || (lo == hi && a <= lo && lo <= b)
                || (k == 0 && hi <= a)
                || (k + 1 == n && lo >= b)
        })
        .collect()
}

/// The curve's `(u, v)` on its wall and its world point as jets in the
/// edge's fraction: on span `pin`'s polynomial, or on the spans the base
/// meets (their union across a knot, up to the order the wall's
/// continuity allows).
pub(super) type WallJet<T> = ([Jet<T>; 2], [Jet<T>; 3]);

pub(super) fn jet<T: Real>(
    m: &WallMeet,
    fraction: &Jet<T>,
    pin: Option<usize>,
) -> Option<WallJet<T>> {
    let sp = spans(m)?;
    let u = fraction.scale(&fc(m.sweep)).add_constant(&fc(m.start));
    let picks = match pin {
        Some(k) => vec![k],
        None => {
            let (lo, hi) = u.c[0].bounds_f64();
            meeting(&sp, lo, hi)
        }
    };
    if picks.is_empty() || (picks.len() > 1 && fraction.order() > sp.continuity + 1) {
        return None;
    }
    let mut out: Option<(Jet<T>, [Jet<T>; 3])> = None;
    for k in picks {
        let (v, p) = eval(m, sp.spans.get(k)?, &u)?;
        out = Some(match out {
            None => (v, p),
            Some((v0, p0)) => (
                union(&v0, &v),
                std::array::from_fn(|i| union(&p0[i], &p[i])),
            ),
        });
    }
    let (v, p) = out?;
    Some(([u, v], p))
}

fn union<T: Real>(a: &Jet<T>, b: &Jet<T>) -> Jet<T> {
    Jet {
        c: a.c.iter().zip(&b.c).map(|(x, y)| x.union(y)).collect(),
    }
}

/// The edge's fraction range split at the wall's knots, exactly: `(fa,
/// fb, span)` with `fa < fb`, ascending.
pub(super) fn pieces(m: &WallMeet) -> Option<Vec<(R, R, usize)>> {
    let sp = spans(m)?;
    let (s, w) = (R::from_float(m.start)?, R::from_float(m.sweep)?);
    if w == zero() {
        return None;
    }
    let e = &s + &w;
    let (a, b) = if w > zero() {
        (s.clone(), e)
    } else {
        (e, s.clone())
    };
    let mut out = Vec::new();
    for (k, span) in sp.spans.iter().enumerate() {
        let lo = if span.u[0] > a {
            span.u[0].clone()
        } else {
            a.clone()
        };
        let hi = if span.u[1] < b {
            span.u[1].clone()
        } else {
            b.clone()
        };
        if lo >= hi {
            continue;
        }
        let (fa, fb) = ((&lo - &s) / &w, (&hi - &s) / &w);
        out.push(if fa < fb { (fa, fb, k) } else { (fb, fa, k) });
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    (!out.is_empty()).then_some(out)
}

/// The edge fractions of the wall's interior knots inside the range
/// (binary64, for pieces of the tessellation's bounds).
pub(crate) fn knot_fractions(m: &WallMeet) -> Vec<f64> {
    let Some(p) = pieces(m) else {
        return Vec::new();
    };
    p.iter()
        .skip(1)
        .map(|(fa, _, _)| crate::solid::split::rational_f64(fa))
        .collect()
}

fn one() -> R {
    R::from_integer(1.into())
}

fn zero() -> R {
    R::from_integer(0.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::certified::Fast;
    use crate::topology::{Curve3, Projection, Surface};
    use crate::{Frame3, KnotVector, Point3, Tolerance, Vec3};

    /// A quadratic wall over two spans (C1 at its knot `u = 1/2`), its
    /// rulings along z over `v` in `[-1, 3]`, met by a cylinder of radius 1.5
    /// about x through `(0, 0, 1)`: its upper or lower curve by `sign`.
    fn wall_meet(sign: f64, start: f64, sweep: f64) -> WallMeet {
        let profile = [(-1.0, 0.0), (-0.3, 0.6), (0.3, 0.4), (1.0, 0.0)];
        let poles = profile
            .iter()
            .flat_map(|&(x, y)| [Point3::new(x, y, -1.0), Point3::new(x, y, 3.0)])
            .collect();
        let wall = crate::BSplineSurface3::new(
            KnotVector::new(2, vec![0.0, 0.5, 1.0], vec![3, 1, 3]).unwrap(),
            KnotVector::new(1, vec![-1.0, 3.0], vec![2, 2]).unwrap(),
            poles,
            None,
        )
        .unwrap();
        let other = Frame3::new(
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Tolerance::default(),
        )
        .unwrap();
        WallMeet {
            wall,
            other,
            other_radius: 1.5,
            sign,
            start,
            sweep,
        }
    }

    fn inside(x: &Fast, v: f64, slack: f64) -> bool {
        let (lo, hi) = x.bounds_f64();
        lo - slack <= v && v <= hi + slack
    }

    /// The binary64 curve lies on the wall and on the cylinder, and its
    /// jets (on the spans a base meets) enclose its parameters and points
    /// and their differences.
    #[test]
    fn meetings_lie_on_the_wall_and_the_cylinder() {
        for sign in [1.0, -1.0] {
            let m = wall_meet(sign, 0.0, 1.0);
            for k in 0..=16 {
                let f = k as f64 / 16.0;
                let p = m.point(f);
                assert!((p.y.hypot(p.z - 1.0) - 1.5).abs() < 1e-14, "{p:?}");
                assert_eq!((p.z - 1.0).signum(), sign);
                let (u, v) = m.parameters(f);
                let q = m.wall.point(u, v).unwrap();
                assert!((p - q).length() < 1e-14, "{p:?} {q:?}");
                let base = Jet::variable(Fast::exact_f64(f), 2);
                let ([ju, jv], jp) = jet(&m, &base, None).unwrap();
                assert!(inside(&ju.c[0], u, 0.0) && inside(&jv.c[0], v, 1e-14));
                let h = 1e-6;
                let (q0, q1) = (m.point(f - h), m.point(f + h));
                for (i, j) in jp.iter().enumerate() {
                    assert!(inside(&j.c[0], p.to_array()[i], 1e-14));
                    let slope = (q1.to_array()[i] - q0.to_array()[i]) / (2.0 * h);
                    let (lo, hi) = j.c[1].bounds_f64();
                    assert!(
                        lo - 1e-5 <= slope && slope <= hi + 1e-5,
                        "{slope} {lo} {hi}"
                    );
                }
            }
        }
    }

    /// Across the knot a jet is the union of both spans' up to the order
    /// the wall's continuity allows (C1: order two), `None` past it; pinned
    /// to one span it is that span's polynomial at any order.
    #[test]
    fn jets_across_a_knot_stop_at_the_walls_continuity() {
        let m = wall_meet(1.0, 0.0, 1.0);
        let at = |n: usize| Jet::variable(Fast::exact_f64(0.5), n);
        assert!(jet(&m, &at(2), None).is_some());
        assert!(jet(&m, &at(3), None).is_none());
        for span in [0, 1] {
            assert!(jet(&m, &at(6), Some(span)).is_some());
        }
        assert!(jet(&m, &Jet::variable(Fast::exact_f64(0.25), 6), None).is_some());
        // The two spans' second derivatives differ at the knot (C1 only).
        let ([_, v0], _) = jet(&m, &at(3), Some(0)).unwrap();
        let ([_, v1], _) = jet(&m, &at(3), Some(1)).unwrap();
        let (a, b) = (v0.c[2].bounds_f64(), v1.c[2].bounds_f64());
        assert!(a.1 < b.0 || b.1 < a.0, "{a:?} {b:?}");
    }

    /// The pieces split the edge's fractions at the knot exactly, ascending
    /// for either direction, their union `[0, 1]`, each on the span holding
    /// its `u`.
    #[test]
    fn pieces_split_at_the_knots_exactly() {
        for (start, sweep) in [(0.1, 0.8), (0.9, -0.8), (0.0, 1.0), (0.6, 0.3)] {
            let m = wall_meet(1.0, start, sweep);
            let p = pieces(&m).unwrap();
            assert_eq!(p[0].0, zero());
            assert_eq!(p[p.len() - 1].1, one());
            for w in p.windows(2) {
                assert_eq!(w[0].1, w[1].0);
            }
            let crosses =
                f64::min(start, start + sweep) < 0.5 && 0.5 < f64::max(start, start + sweep);
            assert_eq!(p.len(), if crosses { 2 } else { 1 });
            let knots = knot_fractions(&m);
            assert_eq!(knots.len(), p.len() - 1);
            for f in knots {
                assert!((start + sweep * f - 0.5).abs() < 1e-15, "{f}");
            }
            let sp = spans(&m).unwrap();
            let (s, w) = (R::from_float(start).unwrap(), R::from_float(sweep).unwrap());
            for (fa, fb, k) in &p {
                let mid = (fa + fb) / R::from_integer(2.into());
                let u = &s + &w * mid;
                assert!(sp.spans[*k].u[0] < u && u < sp.spans[*k].u[1]);
            }
        }
    }

    /// The pcurve on the curve's own wall is the wall's parameters (no
    /// inversion), for either use; on another surface there is none.
    #[test]
    fn the_own_walls_pcurve_is_its_parameters() {
        let m = wall_meet(-1.0, 0.1, 0.8);
        let surface = Surface::BSpline(m.wall.clone());
        for reversed in [false, true] {
            let pr = Projection::own_wall(
                Curve3::WallMeet(Box::new(m.clone())),
                surface.clone(),
                reversed,
            )
            .unwrap();
            for k in 0..=8 {
                let f = k as f64 / 8.0;
                let (u, v) = m.parameters(if reversed { 1.0 - f } else { f });
                let base = Jet::variable(Fast::exact_f64(f), 1);
                let [ju, jv] = super::super::projection::projection_jet(&pr, &base).unwrap();
                assert!(inside(&ju.c[0], u, 1e-15) && inside(&jv.c[0], v, 1e-14));
            }
        }
        let plane = Surface::Plane(Frame3::xy());
        assert!(Projection::own_wall(Curve3::WallMeet(Box::new(m)), plane, false).is_none());
    }
}
