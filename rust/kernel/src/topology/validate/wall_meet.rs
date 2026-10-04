//! Certified evaluation of a spline wall's meeting with a cylinder
//! (`Curve3::WallMeet`, S9f.2b; D13) or a sphere (S9f.3a), in either tier and
//! over any of the integrands' numbers (`quadrature::Num`: plain enclosures,
//! Taylor series, jets).
//!
//! The wall is its face's stored B-spline surface of degree one in `v`: on
//! each knot span its two pole rows are exact Bézier rows (the surface's
//! exact patches), so its ruling there is `L(ū) + t M(ū)` with `L` the low
//! row's Bernstein polynomial, `M` the rows' difference over the `v` range
//! and `ū` the span's local parameter. The other surface's function along
//! the ruling is `a t^2 + 2 b t + c`, and the curve's `t` is `(-b + s
//! sqrt(b^2 - a c)) / a` or `c / (-b - s sqrt(b^2 - a c))`, the one whose
//! middle cancels less (the binary64 curve's own choice, both the same
//! value). Over the other surface's rows `e_i` (a cylinder's two axes, a
//! sphere's three world axes, S9f.3a), with the foot's and the direction's
//! coordinates `w_i` and `m_i` along them, `a = sum m_i^2`, `b = sum w_i
//! m_i`, `c = sum w_i^2 - r^2` and the discriminant `d = a r^2 - sum_{i<j}
//! (w_i m_j - w_j m_i)^2` (Lagrange's identity) are polynomials in `ū` whose
//! Bernstein coefficients are made exactly once per wall and surface
//! (products of Bernstein
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
//!
//! S9f.2b.2: a graph over the wall's `v` (a `WallMeet` with a `window`)
//! lies on one span, its `u` the window's root of the cylinder's function
//! along the ruling (`eval_height`: interval Newton at the base, then the
//! series by Newton's steps about a point and coefficient by coefficient
//! over a range). At a series or a jet every Bernstein polynomial is taken
//! by its Taylor expansion about the base's value (`Coefficients::at`).
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
    /// The Taylor coefficients' Bernstein polynomials, `P^(m) / m! = C(n,
    /// m) sum_i (Δ^m b)_i B_i^(n - m)` for `m >= 1`, exact and enclosed.
    taylor: Vec<(Vec<R>, Vec<Fast>)>,
}

impl Coefficients {
    fn new(exact: Vec<R>) -> Self {
        let fast = exact.iter().map(Fast::from_r).collect();
        let n = exact.len().saturating_sub(1);
        let mut taylor = Vec::new();
        let mut diff = exact.clone();
        let mut binomial = one();
        for m in 1..=n {
            diff = diff.windows(2).map(|w| &w[1] - &w[0]).collect();
            binomial = binomial * R::from_integer(((n - m + 1) as i64).into())
                / R::from_integer((m as i64).into());
            let e: Vec<R> = diff.iter().map(|x| x * &binomial).collect();
            let f = e.iter().map(Fast::from_r).collect();
            taylor.push((e, f));
        }
        Self {
            exact,
            fast,
            taylor,
        }
    }

    /// At an enclosure, by de Casteljau (`bernstein`); at a series or a jet
    /// by its Taylor expansion about the base, `sum_m P^(m)(t_0) / m! (t -
    /// t_0)^m` in Horner's form (`P^(m) / m!` enclosed at the base by de
    /// Casteljau): `n` products of series instead of de Casteljau's `n (n +
    /// 1)` (S9f.2b.2: a graph over `v` evaluates these at every coefficient
    /// it solves for), and over a range a mean-value form, its higher terms
    /// with no constant term to widen.
    fn at<T: Real, N: Num<T>>(&self, t: &N) -> N {
        let scalars = |exact: &[R], fast: &[Fast]| -> Vec<T> {
            exact
                .iter()
                .zip(fast)
                .map(|(x, f)| cached::<T>(f, x))
                .collect()
        };
        let b = scalars(&self.exact, &self.fast);
        if t.terms() <= 1 || self.taylor.is_empty() {
            return bernstein(&b, t);
        }
        let t0 = Jet::constant(t.coefficient_at(0), 0);
        let at0 = |b: &[T]| bernstein(b, &t0).c[0].clone();
        let delta = t.with_coefficient(0, T::exact_f64(0.0));
        let mut acc = t.lift(&at0(&scalars(
            &self.taylor[self.taylor.len() - 1].0,
            &self.taylor[self.taylor.len() - 1].1,
        )));
        for m in (0..self.taylor.len()).rev() {
            let c = if m == 0 {
                at0(&b)
            } else {
                at0(&scalars(&self.taylor[m - 1].0, &self.taylor[m - 1].1))
            };
            acc = acc.mul(&delta).shift(&c);
        }
        acc
    }
}

/// A knot span of the wall met by one surface: its `u` range (exact and
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
    /// The foot's and the direction's coordinates along each of the other
    /// surface's rows (`[w_i, m_i]`), for ranges.
    axes: Vec<[Coefficients; 2]>,
    /// Their derivatives in `ū` (S9f.2b.2's graphs over `v`).
    daxes: Vec<[Coefficients; 2]>,
    /// `-u0`, `1 / (u1 - u0)` and `v0`, exact, and binary64.
    constants: [R; 3],
    fast_constants: [Fast; 3],
}

/// The wall's spans met by the other surface and the wall's least
/// continuity across an interior knot.
pub(super) struct Spans {
    pub(super) spans: Vec<Span>,
    continuity: usize,
}

/// Walls and surfaces kept (each pair's polynomials are exact rational
/// work).
const KEPT: usize = 64;

type Key = (BSplineSurface3, crate::Frame3, u64, bool);

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

/// The Bernstein coefficients of a Bernstein polynomial's derivative (one
/// degree lower; zero for a constant).
fn bernstein_derivative(b: &[R]) -> Vec<R> {
    let n = b.len() - 1;
    if n == 0 {
        return vec![zero()];
    }
    let k = R::from_integer((n as i64).into());
    b.windows(2).map(|w| (&w[1] - &w[0]) * &k).collect()
}

fn combine(f: &[R], g: &[R], scale: &R) -> Vec<R> {
    f.iter().zip(g).map(|(a, b)| a + b * scale).collect()
}

/// The wall's spans met by the curve's other surface, made once per wall
/// and surface (`None` for a wall that is not a nonrational surface of
/// degree one in `v` over one `v` span).
pub(super) fn spans(m: &WallMeet) -> Option<Arc<Spans>> {
    let key: Key = (
        m.wall.clone(),
        m.other,
        m.other_radius.to_bits(),
        m.other_sphere,
    );
    if let Some(hit) = WALLS.with(|w| {
        w.borrow()
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, x)| x.clone())
    }) {
        return Some(hit);
    }
    let exact = |x: f64| R::from_float(x);
    let o = m.other.origin().to_array();
    let o: [R; 3] = [exact(o[0])?, exact(o[1])?, exact(o[2])?];
    // The other surface's rows, exact (a cylinder's stored axes, the
    // world's for a sphere).
    let mut rows: Vec<[R; 3]> = Vec::new();
    for e in m.other_rows() {
        let e = e.to_array();
        rows.push([exact(e[0])?, exact(e[1])?, exact(e[2])?]);
    }
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
        // The foot's and the direction's coordinates along the other
        // surface's rows, then `a`, `b`, `c` and `d`.
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
        let wm: Vec<[Vec<R>; 2]> = rows
            .iter()
            .map(|e| [axis(&low, e, true), axis(&dir, e, false)])
            .collect();
        let minus = -one();
        let zeros = |n: usize| vec![zero(); n];
        let (mut a, mut b, mut c) = (zeros(2 * du + 1), zeros(2 * du + 1), zeros(2 * du + 1));
        for [w, mm] in &wm {
            a = combine(&a, &product(mm, mm), &one());
            b = combine(&b, &product(w, mm), &one());
            c = combine(&c, &product(w, w), &one());
        }
        let c: Vec<R> = c.into_iter().map(|x| x - &r2).collect();
        // `sum_{i<j} (w_i m_j - w_j m_i)^2`, of degree `4 du`.
        let mut cross2 = zeros(4 * du + 1);
        for i in 0..wm.len() {
            for j in i + 1..wm.len() {
                let x = combine(
                    &product(&wm[i][0], &wm[j][1]),
                    &product(&wm[j][0], &wm[i][1]),
                    &minus,
                );
                cross2 = combine(&cross2, &product(&x, &x), &one());
            }
        }
        // `a r^2` raised to the cross terms' degree (times one of degree
        // `2 du`).
        let ar2: Vec<R> = a.iter().map(|x| x * &r2).collect();
        let d = combine(&product(&ar2, &vec![one(); 2 * du + 1]), &cross2, &minus);
        let constants = [-&u[0], one() / (&u[1] - &u[0]), v[0].clone()];
        let fast_constants = std::array::from_fn(|i| Fast::from_r(&constants[i]));
        spans.push(Span {
            uf: [f(&u[0]), f(&u[1])],
            u,
            v,
            low: low.map(Coefficients::new),
            dir: dir.map(Coefficients::new),
            abcd: [a, b, c, d].map(Coefficients::new),
            daxes: wm
                .iter()
                .map(|[w, mm]| [w, mm].map(|x| Coefficients::new(bernstein_derivative(x))))
                .collect(),
            axes: wm.into_iter().map(|x| x.map(Coefficients::new)).collect(),
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
        let wm: Vec<[N; 2]> = span
            .axes
            .iter()
            .map(|[w, mm]| [w.at::<T, N>(&ub), mm.at::<T, N>(&ub)])
            .collect();
        let r = fc::<T>(m.other_radius);
        let r2 = r.mul(&r);
        let a = sum_of(wm.iter().map(|[_, mm]| mm.square()));
        let b = sum_of(wm.iter().map(|[w, mm]| w.mul(mm)));
        let c = sum_of(wm.iter().map(|[w, _]| w.square())).shift(&r2.neg());
        let mut pairs = Vec::new();
        for i in 0..wm.len() {
            for j in i + 1..wm.len() {
                pairs.push(
                    wm[i][0]
                        .mul(&wm[j][1])
                        .sub(&wm[j][0].mul(&wm[i][1]))
                        .square(),
                );
            }
        }
        let d = a.scale(&r2).sub(&sum_of(pairs.into_iter()));
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

/// The sum of a nonempty sequence of numbers.
fn sum_of<T: Real, N: Num<T>>(mut it: impl Iterator<Item = N>) -> N {
    let first = it.next().expect("a row");
    it.fold(first, |acc, x| acc.add(&x))
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
    // S9f.2b.2: a graph over `v` on its window's span.
    if m.window.is_some() {
        let k = window_span(&sp, m)?;
        if pin.is_some_and(|p| p != k) {
            return None;
        }
        let span = &sp.spans[k];
        let v = fraction.scale(&fc(m.sweep)).add_constant(&fc(m.start));
        let t = v.add_constant(&T::from_r(&span.v[0]).neg());
        let (u, _, p) = eval_height(m, span, &t)?;
        return Some(([u, v], p));
    }
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

/// The span holding a graph over `v`'s window (S9f.2b.2).
pub(super) fn window_span(sp: &Spans, m: &WallMeet) -> Option<usize> {
    let [a, b] = m.window?;
    let (a, b) = (R::from_float(a)?, R::from_float(b)?);
    sp.spans
        .iter()
        .position(|s| s.u[0] <= a && b <= s.u[1] && a < b)
}

/// S9f.2b.2: a graph over the wall's `v` at `t = v - v0` (any number of the
/// integrands, its base an enclosure over a range or a point): the root
/// `u` of `g(u, t) = sum X_i^2 - r^2` (`X_i = w_i(ū) + t m_i(ū)` over the
/// other surface's rows: its function along the ruling) in the window, its
/// derivative `du/dt = -g_t / g_u` and the curve's point. The root at the
/// base by interval Newton about the binary64 root (`WallMeet::root_at`),
/// `u* - g(u*, t) / g_u(U, t)` strictly inside `U` and the window (a unique
/// root in `U` at every `t` of the base); its Taylor coefficients by the
/// implicit function theorem term by term: the `k`-th coefficient of
/// `g(u, t)` is `g_u u_k` plus terms of the lower ones, so `u_k` is minus
/// those (the coefficient with `u_k` zero) over `g_u` enclosed over the
/// root's box (inclusion isotone: an enclosure at every point of the base).
pub(super) fn eval_height<T: Real, N: Num<T>>(
    m: &WallMeet,
    span: &Span,
    t: &N,
) -> Option<(N, N, [N; 3])> {
    let window = m.window?;
    let constant = |i: usize| cached::<T>(&span.fast_constants[i], &span.constants[i]);
    let local = |u: &N| u.shift(&constant(0)).scale(&constant(1));
    let r2 = fc::<T>(m.other_radius).square();
    // Each row's coordinate `X_i = w_i + t m_i` and its derivative in `u`
    // at `u` (any number).
    let parts = |u: &N, t: &N| -> Vec<[N; 2]> {
        let ub = local(u);
        span.axes
            .iter()
            .zip(&span.daxes)
            .map(|([w, mm], [dw, dm])| {
                let x = w.at::<T, N>(&ub).add(&t.mul(&mm.at::<T, N>(&ub)));
                let xu = dw
                    .at::<T, N>(&ub)
                    .add(&t.mul(&dm.at::<T, N>(&ub)))
                    .scale(&constant(1));
                [x, xu]
            })
            .collect()
    };
    // `g` and `g_u` over enclosures (a jet of order one in `u`).
    let g_gu = |u: &T, t: &T| -> (T, T) {
        let uj = Jet::variable(u.clone(), 1);
        let tj = Jet::constant(t.clone(), 1);
        let ub = uj.add_constant(&constant(0)).scale(&constant(1));
        let g = sum_of(span.axes.iter().map(|[w, mm]| {
            w.at::<T, Jet<T>>(&ub)
                .add(&tj.mul(&mm.at::<T, Jet<T>>(&ub)))
                .square()
        }))
        .shift(&r2.neg());
        (g.c[0].clone(), g.c[1].clone())
    };
    let tb = t.coefficient_at(0);
    let (tlo, thi) = tb.bounds_f64();
    let star = m.root_at(0.5 * tlo + 0.5 * thi);
    if !star.is_finite() || star <= window[0] || star >= window[1] {
        return None;
    }
    let s0 = fc::<T>(star);
    let (g0, _) = g_gu(&s0, &tb);
    let least = |x: &T| {
        let (a, b) = x.bounds_f64();
        if a > 0.0 {
            a
        } else if b < 0.0 {
            -b
        } else {
            0.0
        }
    };
    let size = |x: &T| {
        let (a, b) = x.bounds_f64();
        a.abs().max(b.abs())
    };
    let (_, d0) = g_gu(&s0, &tb);
    let slope0 = least(&d0);
    if slope0 <= 0.0 || slope0.is_nan() {
        return None;
    }
    let mut delta = 4.0 * size(&g0) / slope0 + 1e-15 * (1.0 + star.abs());
    let mut root = None;
    for _ in 0..12 {
        if !delta.is_finite() || star - delta <= window[0] || star + delta >= window[1] {
            break;
        }
        let box_u = fc::<T>(star - delta).union(&fc(star + delta));
        let (_, d) = g_gu(&box_u, &tb);
        if let Some(q) = g0.div(&d) {
            let next = s0.sub(&q);
            let (nlo, nhi) = next.bounds_f64();
            if nlo > star - delta && nhi < star + delta {
                root = Some(next);
                break;
            }
        }
        delta *= 4.0;
    }
    let u0 = root?;
    let (_, gu0) = g_gu(&u0, &tb);
    if least(&gu0) <= 0.0 {
        return None;
    }
    // About a point, Newton's steps on the series, each doubling the
    // coefficients known: with `u`'s first `k` those of the root (enclosed)
    // and the rest any, `u - g(u, t) / g_u(u, t)` holds the root's first `2
    // k`, its constant term kept the interval Newton enclosure (the
    // recurrence below divides by `g_u` once per coefficient, widening a
    // point's high coefficients by its rounding each time). Over a range,
    // coefficient by coefficient: the `k`-th of `g(u, t)` with `u_k` zero
    // over `g_u` enclosed over the root's box (Newton's quotient of two wide
    // series overestimates by orders of magnitude there). Both inclusion
    // isotone: enclosed at every point of the base.
    let mut u = t.lift(&u0);
    if t.sharp() {
        let mut known = 1;
        while known < t.terms() {
            let rows = parts(&u, t);
            let g = sum_of(rows.iter().map(|[x, _]| x.square())).shift(&r2.neg());
            let gu = sum_of(rows.iter().map(|[x, xu]| x.mul(xu))).scale(&fc(2.0));
            let next = u.sub(&g.div(&gu)?);
            known = (2 * known).min(t.terms());
            for k in 1..known {
                u = u.with_coefficient(k, next.coefficient_at(k));
            }
        }
    } else {
        for k in 1..t.terms() {
            let ub = local(&u);
            let g = sum_of(
                span.axes
                    .iter()
                    .map(|[w, mm]| w.at::<T, N>(&ub).add(&t.mul(&mm.at::<T, N>(&ub))).square()),
            )
            .shift(&r2.neg());
            let uk = g.coefficient_at(k).div(&gu0)?.neg();
            u = u.with_coefficient(k, uk);
        }
    }
    let ub = local(&u);
    let rows = parts(&u, t);
    let gt = sum_of(
        rows.iter()
            .zip(&span.axes)
            .map(|([x, _], [_, mm])| x.mul(&mm.at::<T, N>(&ub))),
    );
    let gu = sum_of(rows.iter().map(|[x, xu]| x.mul(xu)));
    let du = gt.div(&gu)?.neg();
    let low: [N; 3] = std::array::from_fn(|k| span.low[k].at(&ub));
    let dir: [N; 3] = std::array::from_fn(|k| span.dir[k].at(&ub));
    let point = std::array::from_fn(|k| low[k].add(&dir[k].mul(t)));
    Some((u, du, point))
}

/// The edge's fraction range split at the wall's knots, exactly: `(fa,
/// fb, span)` with `fa < fb`, ascending; a graph over `v` (S9f.2b.2) one
/// piece on its window's span.
pub(super) fn pieces(m: &WallMeet) -> Option<Vec<(R, R, usize)>> {
    let sp = spans(m)?;
    if m.window.is_some() {
        return Some(vec![(zero(), one(), window_span(&sp, m)?)]);
    }
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
            other_sphere: false,
            sign,
            start,
            sweep,
            window: None,
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

    /// S9f.2b.2: the wall above met by a thin rod about `x` through `(0,
    /// 0.3, 1)` of radius 0.2, whose meeting turns back where the profile's
    /// `y` is 0.1 (the ruling tangent to the rod): a graph over `v` about
    /// that turning point, its window across it. Its binary64 points lie on
    /// the wall and the rod at the window's root, its jets (over a point
    /// and over a range) enclose its parameters, points and slopes, and its
    /// pieces are one, on the window's span.
    #[test]
    fn graphs_over_v_lie_on_the_wall_and_the_cylinder() {
        let base = wall_meet(1.0, 0.0, 1.0);
        let other = Frame3::new(
            Point3::new(0.0, 0.3, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Tolerance::default(),
        )
        .unwrap();
        let y = |u: f64| base.wall.point(u, 0.0).unwrap().y;
        // The turning point on the first span: y = 0.1 rising.
        let (mut a, mut b) = (0.0, 0.25);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            if y(mid) < 0.1 {
                a = mid;
            } else {
                b = mid;
            }
        }
        let turn = 0.5 * (a + b);
        assert!(turn > 0.01 && turn < 0.2, "{turn}");
        let switch = turn + 0.01;
        let h = (0.04 - (y(switch) - 0.3).powi(2)).sqrt();
        // `v` is the height `z` here: the rod's axis at v = 1.
        let m = WallMeet {
            other,
            other_radius: 0.2,
            window: Some([turn - 0.01, turn + 0.015]),
            start: 1.0 - h,
            sweep: 2.0 * h,
            ..base
        };
        for k in 0..=16 {
            let f = k as f64 / 16.0;
            let p = m.point(f);
            assert!(((p.y - 0.3).hypot(p.z - 1.0) - 0.2).abs() < 1e-14, "{p:?}");
            let (u, v) = m.parameters(f);
            assert!(u > turn - 0.01 && u < turn + 0.015);
            let q = m.wall.point(u, v).unwrap();
            assert!((p - q).length() < 1e-14, "{p:?} {q:?}");
            for at in [
                Jet::variable(Fast::exact_f64(f), 3),
                Jet::variable(Fast::exact_f64(f).union(&Fast::exact_f64(f + 1e-3)), 3),
            ] {
                let ([ju, jv], jp) = jet(&m, &at, None).unwrap();
                assert!(inside(&ju.c[0], u, 1e-15) && inside(&jv.c[0], v, 1e-14));
                let dh = 1e-6;
                let (q0, q1) = (m.point(f - dh), m.point(f + dh));
                for (i, j) in jp.iter().enumerate() {
                    assert!(inside(&j.c[0], p.to_array()[i], 1e-14));
                    let slope = (q1.to_array()[i] - q0.to_array()[i]) / (2.0 * dh);
                    let (lo, hi) = j.c[1].bounds_f64();
                    assert!(
                        lo - 1e-5 <= slope && slope <= hi + 1e-5,
                        "{slope} {lo} {hi}"
                    );
                }
            }
        }
        // Over a range the coefficients hold the point's and stay narrow
        // (term by term: Newton's series quotient over this range left the
        // thirteenth `±0.3`, here `±1.6e-3`); about a point the high ones
        // stay thin (Newton's steps: term by term left the thirteenth `±1.2`).
        let base = Fast::exact_f64(0.475).union(&Fast::exact_f64(0.525));
        let range = jet(&m, &Jet::variable(base, 13), None).unwrap().1[1].clone();
        let point = jet(&m, &Jet::variable(Fast::exact_f64(0.5), 13), None)
            .unwrap()
            .1[1]
            .clone();
        for k in 0..=13 {
            let (a, b) = point.c[k].bounds_f64();
            let (c, d) = range.c[k].bounds_f64();
            assert!(c <= a && b <= d, "{k}: [{a}, {b}] in [{c}, {d}]");
            assert!(d - c < 1e-2 && b - a < 1e-9, "{k}: [{c}, {d}], [{a}, {b}]");
        }
        assert_eq!(pieces(&m).unwrap(), vec![(zero(), one(), 0)]);
        assert!(knot_fractions(&m).is_empty());
        // A window across the knot holds no span.
        let across = WallMeet {
            window: Some([0.4, 0.6]),
            ..m.clone()
        };
        assert!(pieces(&across).is_none());
    }

    /// Jets of order `n` over a point and over a short range enclose the
    /// binary64 curve's parameters, points and slopes at `f` (across a knot
    /// of the C1 wall order two at most).
    fn jets_enclose(m: &WallMeet, f: f64, n: usize, name: &str) {
        let (u, v) = m.parameters(f);
        let p = m.point(f);
        for at in [
            Jet::variable(Fast::exact_f64(f), n),
            Jet::variable(Fast::exact_f64(f).union(&Fast::exact_f64(f + 1e-3)), n),
        ] {
            let ([ju, jv], jp) = jet(m, &at, None).unwrap();
            assert!(
                inside(&ju.c[0], u, 1e-14) && inside(&jv.c[0], v, 1e-14),
                "{name}"
            );
            let dh = 1e-6;
            let (q0, q1) = (m.point(f - dh), m.point(f + dh));
            for (i, j) in jp.iter().enumerate() {
                assert!(inside(&j.c[0], p.to_array()[i], 1e-14), "{name}");
                let slope = (q1.to_array()[i] - q0.to_array()[i]) / (2.0 * dh);
                let (lo, hi) = j.c[1].bounds_f64();
                assert!(
                    lo - 1e-5 <= slope && slope <= hi + 1e-5,
                    "{name}: {slope} {lo} {hi}"
                );
            }
        }
    }

    /// S9f.3a: the wall above met by a sphere about `(0.1, 0.3, 1)` of
    /// radius 1.2 (every ruling twice: the profile within 1.14 of the
    /// centre's foot), its upper or lower curve by `sign` over the whole
    /// `u` domain; and a small sphere about `(0, 0.4, 1)` of radius 0.2
    /// turning back where the profile is 0.2 from its centre's foot, a
    /// graph over `v` about that turning point. Their binary64 points lie on
    /// the wall and the sphere, their jets enclose them (the sphere's three
    /// rows in `a`, `b`, `c` and `d`).
    #[test]
    fn sphere_meetings_lie_on_the_wall_and_the_sphere() {
        let at = |x: f64, y: f64, z: f64| {
            Frame3::new(
                Point3::new(x, y, z),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                Tolerance::default(),
            )
            .unwrap()
        };
        for sign in [1.0, -1.0] {
            let m = WallMeet {
                other: at(0.1, 0.3, 1.0),
                other_radius: 1.2,
                other_sphere: true,
                ..wall_meet(sign, 0.0, 1.0)
            };
            let c = m.other.origin();
            for k in 0..=16 {
                let f = f64::from(k) / 16.0;
                let p = m.point(f);
                assert!(((p - c).length() - 1.2).abs() < 1e-14, "{p:?}");
                assert_eq!((p.z - 1.0).signum(), sign);
                let (u, v) = m.parameters(f);
                assert!((m.wall.point(u, v).unwrap() - p).length() < 1e-14);
                if k < 16 {
                    jets_enclose(&m, f, 2, "over u");
                }
            }
        }
        let base = wall_meet(1.0, 0.0, 1.0);
        let reach = |u: f64| {
            let p = base.wall.point(u, 0.0).unwrap();
            p.x.hypot(p.y - 0.4)
        };
        // The turning point on the first span: the profile 0.2 from (0, 0.4).
        let (mut a, mut b) = (0.0, 0.5);
        assert!(reach(a) > 0.2 && reach(b) < 0.2);
        for _ in 0..100 {
            let mid = 0.5 * (a + b);
            if reach(mid) > 0.2 {
                a = mid;
            } else {
                b = mid;
            }
        }
        let turn = 0.5 * (a + b);
        let switch = turn + 0.01;
        let h = (0.04 - reach(switch).powi(2)).sqrt();
        let m = WallMeet {
            other: at(0.0, 0.4, 1.0),
            other_radius: 0.2,
            other_sphere: true,
            window: Some([turn - 0.01, turn + 0.015]),
            start: 1.0 - h,
            sweep: 2.0 * h,
            ..base
        };
        let c = m.other.origin();
        for k in 0..16 {
            let f = f64::from(k) / 16.0;
            let p = m.point(f);
            assert!(((p - c).length() - 0.2).abs() < 1e-14, "{p:?}");
            let (u, v) = m.parameters(f);
            assert!(u > turn - 0.01 && u < turn + 0.015);
            assert!((m.wall.point(u, v).unwrap() - p).length() < 1e-14);
            jets_enclose(&m, f, 3, "over v");
        }
    }
}
