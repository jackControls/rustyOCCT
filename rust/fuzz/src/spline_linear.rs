//! Complete spline/line and spline/segment preimages from constructed answers.
//! Known factors, per-span rational linear inequalities on polylines, and
//! closed-form quadratic sign functions define every expected parameter and
//! interval. Production instead isolates common gcd roots and clips by Taylor
//! signs at algebraic boundaries; no expected value comes from its isolator.
use crate::byte;
use crate::splines::{bounds, rat};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::*;
use rusty_occt::{Error, ExactBSplineCurve3, ExactKnotVector, ScalarInterval};
use std::cmp::Ordering::{self, Equal, Greater, Less};
use std::rc::Rc;

fn r(n: i64) -> R {
    R::from_integer(n.into())
}
fn q(n: i64, d: i64) -> R {
    R::new(n.into(), d.into())
}
fn pow2(n: i32) -> R {
    let x = BigInt::from(1) << n.unsigned_abs() as usize;
    if n < 0 {
        R::new(1.into(), x)
    } else {
        R::from_integer(x)
    }
}
fn choose(n: usize, k: usize) -> R {
    R::from_integer((0..k).fold(BigInt::from(1), |a, i| a * (n - i) / (i + 1)))
}
fn on_axis(x: R) -> [R; 3] {
    [x, r(0), r(0)]
}

struct Bytes<'a>(&'a [u8], usize);
impl Bytes<'_> {
    fn next(&mut self) -> u8 {
        self.1 += 1;
        byte(self.0, self.1 - 1)
    }
    /// Uniform-ish integer in -offset..n-offset.
    fn small(&mut self, n: u8, offset: i64) -> i64 {
        i64::from(self.next() % n) - offset
    }
}

/// A real number known exactly, or by its exact sign against every rational
/// probe. `cmp(v)` is value.cmp(v).
#[derive(Clone)]
enum Value {
    Exact(R),
    Algebraic(Rc<dyn Fn(&R) -> Ordering>),
}
impl Value {
    fn cmp(&self, v: &R) -> Ordering {
        match self {
            Value::Exact(x) => x.cmp(v),
            Value::Algebraic(f) => f(v),
        }
    }
    /// offset + scale*self, for any rational scale.
    fn affine(&self, offset: &R, scale: &R) -> Value {
        match self {
            Value::Exact(x) => Value::Exact(offset + scale * x),
            _ if scale == &r(0) => Value::Exact(offset.clone()),
            Value::Algebraic(f) => {
                let (f, offset, scale) = (f.clone(), offset.clone(), scale.clone());
                Value::Algebraic(Rc::new(move |v| {
                    let order = f(&((v - &offset) / &scale));
                    if scale > r(0) {
                        order
                    } else {
                        order.reverse()
                    }
                }))
            }
        }
    }
}

/// sqrt(g).cmp(v) for rational g >= 0.
fn sqrt_cmp(g: &R, v: &R) -> Ordering {
    if v < &r(0) {
        Greater
    } else {
        g.cmp(&(v * v))
    }
}
fn rational_sqrt(g: &R) -> Option<R> {
    let (n, d) = (g.numer().sqrt(), g.denom().sqrt());
    (&n * &n == *g.numer() && &d * &d == *g.denom()).then(|| R::new(n, d))
}

/// Curve parameter t (before the a+width*t domain map), base point
/// `scale*direction`, and linear query parameter.
#[derive(Clone)]
struct Hit {
    t: Value,
    scale: Value,
    direction: [R; 3],
    linear: Value,
}
fn rational_hit(t: R, point: [R; 3], linear: R) -> Hit {
    Hit {
        t: Value::Exact(t),
        scale: Value::Exact(r(1)),
        direction: point,
        linear: Value::Exact(linear),
    }
}

struct Case {
    degree: usize,
    periodic: bool,
    /// Homogeneous controls in base coordinates.
    controls: Vec<[R; 4]>,
    /// Knots and closed query range in t units.
    knots: Vec<R>,
    mults: Vec<usize>,
    lo: R,
    hi: R,
    /// Base query endpoints.
    a: [R; 3],
    b: [R; 3],
    points: Vec<Hit>,
    intervals: Vec<[Hit; 2]>,
}

/// Base-axis endpoints; segments may collapse or reverse, lines may not collapse.
fn axis_endpoints(bytes: &mut Bytes, line: bool) -> (R, R) {
    let a = q(bytes.small(13, 6), 4);
    let mut b = if bytes.next() % 4 == 0 {
        a.clone()
    } else {
        q(bytes.small(13, 6), 4)
    };
    if line && a == b {
        b = &a + r(1);
    }
    (a, b)
}
/// Linear parameter of x on the axis query, when it is in the query set.
fn axis_parameter(x: &R, a: &R, b: &R, line: bool) -> Option<R> {
    let s = if a == b {
        (x == a).then(|| r(0))?
    } else {
        (x - a) / (b - a)
    };
    (line || (r(0) <= s && s <= r(1))).then_some(s)
}

/// Degree n Bernstein curve (t/W, 0, N/W) where N has prescribed rational roots.
fn known_factor(bytes: &mut Bytes, line: bool) -> Case {
    let degree = 1 + usize::from(bytes.next() % 25);
    let denominator = 2 + bytes.small(8, 0);
    let mut roots: Vec<R> = (0..degree)
        .map(|_| q(bytes.small((denominator + 3) as u8, 1), denominator))
        .collect();
    if degree >= 2 && bytes.next() % 4 == 0 {
        roots[0] = q(1, 3);
        roots[1] = q(1, 3) + pow2(-(1 + i32::from(bytes.next())) * 4);
    }
    let weights: Vec<R> = (0..=degree)
        .map(|_| q(1 + bytes.small(31, 0), 11))
        .collect();
    let weight = |t: &R| -> R {
        (0..=degree)
            .map(|i| {
                &weights[i]
                    * choose(degree, i)
                    * t.pow(i as i32)
                    * (r(1) - t).pow((degree - i) as i32)
            })
            .sum()
    };
    let mut polynomial = vec![r(1)];
    for root in &roots {
        let mut next = vec![r(0); polynomial.len() + 1];
        for (i, c) in polynomial.iter().enumerate() {
            next[i] -= c * root;
            next[i + 1] += c;
        }
        polynomial = next;
    }
    let controls = (0..=degree)
        .map(|i| {
            let n: R = (0..=i)
                .map(|j| &polynomial[j] * choose(i, j) / choose(degree, j))
                .sum();
            [R::new(i.into(), degree.into()), r(0), n, weights[i].clone()]
        })
        .collect();
    let lo = q(bytes.small(9, 0), 8);
    let hi = &lo + (r(1) - &lo) * q(bytes.small(9, 0), 8);
    let mut distinct = roots;
    distinct.sort();
    distinct.dedup();
    distinct.retain(|t| lo <= *t && *t <= hi);
    let (mut a, mut b) = axis_endpoints(bytes, line);
    if !line && bytes.next() % 4 == 0 {
        // Collapse the segment onto one known contact.
        if let Some(t) = distinct.first() {
            a = t / weight(t);
            b = a.clone();
        }
    }
    let points = distinct
        .iter()
        .filter_map(|t| {
            let x = t / weight(t);
            let s = axis_parameter(&x, &a, &b, line)?;
            Some(rational_hit(t.clone(), on_axis(x), s))
        })
        .collect();
    Case {
        degree,
        periodic: false,
        controls,
        knots: vec![r(0), r(1)],
        mults: vec![degree + 1; 2],
        lo,
        hi,
        a: on_axis(a),
        b: on_axis(b),
        points,
        intervals: vec![],
    }
}

/// Degree-one rational polyline, optionally periodic with several turns.
/// Every span is linear in s, so every expected endpoint is rational.
fn polyline(bytes: &mut Bytes, line: bool) -> Case {
    let periodic = bytes.next() % 3 == 0;
    let spans = 1 + usize::from(bytes.next() % 6) + usize::from(periodic);
    let vertices = spans + usize::from(!periodic);
    let off_axis = bytes.next() % 5;
    let mut controls: Vec<[R; 4]> = (0..vertices)
        .map(|_| {
            let x = q(bytes.small(9, 4), 2);
            let (mut y, mut z) = (r(0), r(0));
            match bytes.next() % 8 {
                v if v < off_axis => y = r(bytes.small(5, 2)),
                7 => z = r(bytes.small(3, 1)),
                _ => {}
            }
            let w = q(1 + bytes.small(7, 0), 3);
            [x * &w, y * &w, z * &w, w]
        })
        .collect();
    if bytes.next() % 4 == 0 {
        // A repeated Cartesian vertex with another weight gives a constant span.
        let i = usize::from(bytes.next()) % (vertices - 1);
        let factor = q(1 + bytes.small(3, 0), 2);
        controls[i + 1] = controls[i].clone().map(|c| c * &factor);
    }
    let n = r(spans as i64);
    let (mults, lo, hi) = if periodic {
        let lo = q(bytes.small(25, 12), 16) * &n;
        let hi = &lo + q(bytes.small(33, 0), 8) * &n;
        (vec![1; spans + 1], lo, hi)
    } else {
        let lo = q(bytes.small(9, 0), 8) * &n;
        let hi = &lo + (&n - &lo) * q(bytes.small(9, 0), 8);
        let mut mults = vec![1; spans + 1];
        mults[0] = 2;
        mults[spans] = 2;
        (mults, lo, hi)
    };
    let (a, b) = axis_endpoints(bytes, line);
    let (low, high) = if a <= b { (&a, &b) } else { (&b, &a) };
    let turns: Vec<i64> = if periodic {
        let turn = |x: &R| -> i64 { (x / &n).floor().to_integer().try_into().unwrap() };
        (turn(&lo) - 1..=turn(&hi) + 1).collect()
    } else {
        vec![0]
    };
    let mut pieces: Vec<[R; 2]> = vec![];
    for turn in turns {
        for j in 0..spans {
            let (p0, p1) = (&controls[j], &controls[(j + 1) % vertices]);
            // Homogeneous component k is c0 + c1 s for s in [0,1].
            let linear = |k: usize| (p0[k].clone(), &p1[k] - &p0[k]);
            let mut range = [r(0), r(1)];
            let mut feasible = true;
            // Y = 0 and Z = 0 on the x axis.
            for (c0, c1) in [linear(1), linear(2)] {
                if c1 == r(0) {
                    feasible &= c0 == r(0);
                } else {
                    let s = -c0 / c1;
                    range[0] = range[0].clone().max(s.clone());
                    range[1] = range[1].clone().min(s);
                }
            }
            // low*W <= X <= high*W with W > 0.
            if !line {
                let (x, w) = (linear(0), linear(3));
                for (c0, c1) in [
                    (&x.0 - low * &w.0, &x.1 - low * &w.1),
                    (high * &w.0 - &x.0, high * &w.1 - &x.1),
                ] {
                    match c1.cmp(&r(0)) {
                        Equal => feasible &= c0 >= r(0),
                        Greater => range[0] = range[0].clone().max(-&c0 / &c1),
                        Less => range[1] = range[1].clone().min(-&c0 / &c1),
                    }
                }
            }
            let base = r(turn * spans as i64 + j as i64);
            let t0 = (&base + &range[0]).max(lo.clone());
            let t1 = (&base + &range[1]).min(hi.clone());
            if feasible && range[0] <= range[1] && t0 <= t1 {
                pieces.push([t0, t1]);
            }
        }
    }
    pieces.sort();
    let mut merged: Vec<[R; 2]> = vec![];
    for piece in pieces {
        match merged.last_mut() {
            Some(last) if last[1] >= piece[0] => {
                if last[1] < piece[1] {
                    last[1] = piece[1].clone();
                }
            }
            _ => merged.push(piece),
        }
    }
    let hit = |t: &R| {
        let local = if periodic {
            t - (t / &n).floor() * &n
        } else {
            t.clone()
        };
        let j = usize::try_from(local.floor().to_integer())
            .unwrap()
            .min(spans - 1);
        let s = &local - r(j as i64);
        let (p0, p1) = (&controls[j], &controls[(j + 1) % vertices]);
        let h: [R; 4] = std::array::from_fn(|k| &p0[k] + &s * (&p1[k] - &p0[k]));
        let point: [R; 3] = std::array::from_fn(|k| &h[k] / &h[3]);
        let linear = axis_parameter(&point[0], &a, &b, true).unwrap();
        rational_hit(t.clone(), point, linear)
    };
    let (mut points, mut intervals) = (vec![], vec![]);
    for [t0, t1] in merged {
        if t0 == t1 {
            points.push(hit(&t0));
        } else {
            intervals.push([hit(&t0), hit(&t1)]);
        }
    }
    Case {
        degree: 1,
        periodic,
        controls,
        knots: (0..=spans).map(|i| r(i as i64)).collect(),
        mults,
        lo,
        hi,
        a: on_axis(a),
        b: on_axis(b),
        points,
        intervals,
    }
}

/// Rational quarter circle (1-t^2, 2t, 0)/(1+t^2) against y = m x through 0.
/// The contact t* is the nonnegative root of m t^2 + 2t - m.
fn circle(bytes: &mut Bytes, line: bool) -> Case {
    let m = q(bytes.small(15, 0), 1 + bytes.small(7, 0));
    let c = match q(bytes.small(17, 4), 1 + bytes.small(8, 0)) {
        c if c == r(0) => r(1),
        c => c,
    };
    let lo = q(bytes.small(9, 0), 16);
    let hi = &lo + (r(1) - &lo) * q(bytes.small(9, 0), 8);
    let one_plus = r(1) + &m * &m;
    // x* = 1/sqrt(1+m^2): x* > v iff v <= 0 or v^2 (1+m^2) < 1.
    let x = {
        let one_plus = one_plus.clone();
        Value::Algebraic(Rc::new(move |v: &R| {
            if v <= &r(0) {
                Greater
            } else {
                r(1).cmp(&(v * v * &one_plus))
            }
        }))
    };
    let t = if m == r(0) {
        Value::Exact(r(0))
    } else {
        // m t^2 + 2t - m increases on t >= 0 and is negative at 0.
        let m = m.clone();
        Value::Algebraic(Rc::new(move |p: &R| {
            if p < &r(0) {
                Greater
            } else {
                r(0).cmp(&(&m * p * p + r(2) * p - &m))
            }
        }))
    };
    let accepted = t.cmp(&lo) != Less
        && t.cmp(&hi) != Greater
        && (line || (c > r(0) && &c * &c * &one_plus >= r(1)));
    let points = if accepted {
        vec![Hit {
            t,
            scale: x.clone(),
            direction: [r(1), m.clone(), r(0)],
            linear: x.affine(&r(0), &(r(1) / &c)),
        }]
    } else {
        vec![]
    };
    Case {
        degree: 2,
        periodic: false,
        controls: vec![
            [r(1), r(0), r(0), r(1)],
            [r(1), r(1), r(0), r(1)],
            [r(0), r(2), r(0), r(2)],
        ],
        knots: vec![r(0), r(1)],
        mults: vec![3, 3],
        lo,
        hi,
        a: on_axis(r(0)),
        b: [c.clone(), &m * &c, r(0)],
        points,
        intervals: vec![],
    }
}

/// x(t) = (2t-1)^2 retraces the x axis. Values g in the query have parameters
/// (1 +- sqrt g)/2, generally irrational interval endpoints.
fn retrace(bytes: &mut Bytes, line: bool) -> Case {
    let alpha = q(bytes.small(6, 0), 4);
    let beta = if !line && bytes.next() % 4 == 0 {
        alpha.clone()
    } else {
        &alpha + q(1 + bytes.small(5, 0), 4)
    };
    let (a, b) = if bytes.next() & 1 == 0 {
        (alpha.clone(), beta.clone())
    } else {
        (beta.clone(), alpha.clone())
    };
    let endpoint = |g: &R, sign: i64| -> Hit {
        let t = match rational_sqrt(g) {
            Some(root) => Value::Exact((r(1) + r(sign) * root) / r(2)),
            None => {
                let g = g.clone();
                // (1 + sign*sqrt g)/2 vs p  <=>  sign*sqrt g vs 2p-1.
                Value::Algebraic(Rc::new(move |p: &R| {
                    let d = r(2) * p - r(1);
                    if sign > 0 {
                        sqrt_cmp(&g, &d)
                    } else {
                        sqrt_cmp(&g, &-d).reverse()
                    }
                }))
            }
        };
        let linear = axis_parameter(g, &a, &b, true).unwrap();
        Hit {
            t,
            scale: Value::Exact(g.clone()),
            direction: on_axis(r(1)),
            linear: Value::Exact(linear),
        }
    };
    let (mut points, mut intervals) = (vec![], vec![]);
    let mut add = |x: Hit, y: Hit| {
        if let (Value::Exact(u), Value::Exact(v)) = (&x.t, &y.t) {
            if u == v {
                points.push(x);
                return;
            }
        }
        intervals.push([x, y]);
    };
    if line {
        add(endpoint(&r(1), -1), endpoint(&r(1), 1));
    } else if alpha <= r(1) {
        let top = beta.clone().min(r(1));
        if alpha == r(0) {
            add(endpoint(&top, -1), endpoint(&top, 1));
        } else if alpha == top {
            points.push(endpoint(&alpha, -1));
            points.push(endpoint(&alpha, 1));
        } else {
            add(endpoint(&top, -1), endpoint(&alpha, -1));
            add(endpoint(&alpha, 1), endpoint(&top, 1));
        }
    }
    Case {
        degree: 2,
        periodic: false,
        controls: [1, -1, 1].map(|x| [r(x), r(0), r(0), r(1)]).to_vec(),
        knots: vec![r(0), r(1)],
        mults: vec![3, 3],
        lo: r(0),
        hi: r(1),
        a: on_axis(a),
        b: on_axis(b),
        points,
        intervals,
    }
}

/// Minimal finite enclosure of an exact or algebraic expected value, or a
/// typed conversion failure outside binary64 range.
fn enclosure(actual: rusty_occt::Result<ScalarInterval>, expected: &Value) {
    let max = rat(f64::MAX);
    if expected.cmp(&max) == Greater || expected.cmp(&-&max) == Less {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
        return;
    }
    let b = actual.unwrap();
    match expected {
        Value::Exact(x) => bounds(b, x),
        Value::Algebraic(_) => {
            let (lo, hi) = (b.lower(), b.upper());
            if lo == hi {
                assert_eq!(expected.cmp(&rat(lo)), Equal);
            } else {
                assert_eq!(lo.next_up(), hi);
                assert_eq!(expected.cmp(&rat(lo)), Greater);
                assert_eq!(expected.cmp(&rat(hi)), Less);
            }
        }
    }
}
fn representable(expected: &Value) -> bool {
    let max = rat(f64::MAX);
    expected.cmp(&max) != Greater && expected.cmp(&-&max) != Less
}

/// Check the kernel's exact comparison against offset + scale*base. Algebraic
/// values are bracketed to 2^-96 in base units by independent bisection.
fn identity(what: &str, base: &Value, offset: &R, scale: &R, compare: impl Fn(&R) -> Ordering) {
    let map = |x: &R| offset + scale * x;
    let mut lo = -r(1);
    let mut hi = r(1);
    if let Value::Exact(x) = base {
        assert_eq!(compare(&map(x)), Equal, "{what}");
        return;
    }
    if scale == &r(0) {
        assert_eq!(compare(offset), Equal, "{what}");
        return;
    }
    while base.cmp(&hi) != Less {
        hi *= r(2);
    }
    while base.cmp(&lo) != Greater {
        lo *= r(2);
    }
    for _ in 0..96 + 8 {
        let mid = (&lo + &hi) / r(2);
        match base.cmp(&mid) {
            Greater => lo = mid,
            Less => hi = mid,
            // An algebraic description can still denote a rational.
            Equal => {
                assert_eq!(compare(&map(&mid)), Equal, "{what}");
                return;
            }
        }
    }
    let (below, above) = if scale > &r(0) {
        (Greater, Less)
    } else {
        (Less, Greater)
    };
    assert_eq!(compare(&map(&lo)), below, "{what}: {lo}");
    assert_eq!(compare(&map(&hi)), above, "{what}: {hi}");
}

pub fn check_spline_linear(data: &[u8]) {
    let mut bytes = Bytes(data, 0);
    let family = bytes.next() % 4;
    let line = bytes.next() % 3 == 0;
    let case = match family {
        0 => known_factor(&mut bytes, line),
        1 => polyline(&mut bytes, line),
        2 => circle(&mut bytes, line),
        _ => retrace(&mut bytes, line),
    };
    // Parameter domain map u = offset + width*t, including far outside binary64.
    let (offset, width) = match bytes.next() % 7 {
        0 => (q(1, 3), q(5, 7)),
        1 => (pow2(2048), pow2(-2048)),
        2 => (-pow2(2048), q(7, 3)),
        3 => (-pow2(-2048), pow2(-2047)),
        4 => (
            q(bytes.small(255, 127), 257),
            pow2(i32::from(bytes.next() as i8)),
        ),
        5 => (pow2(1024), pow2(1024)),
        _ => (r(0), r(1)),
    };
    let global = |t: &R| &offset + &width * t;
    // An invertible integer affine map preserves incidence and the linear
    // query parameter; it turns the axis-aligned constructions oblique.
    let mut matrix: [[R; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| r(i64::from(i == j))));
    if bytes.next() & 1 == 1 {
        let candidate: [[R; 3]; 3] =
            std::array::from_fn(|_| std::array::from_fn(|_| r(bytes.small(7, 3))));
        let m = &candidate;
        let det = &m[0][0] * (&m[1][1] * &m[2][2] - &m[1][2] * &m[2][1])
            - &m[0][1] * (&m[1][0] * &m[2][2] - &m[1][2] * &m[2][0])
            + &m[0][2] * (&m[1][0] * &m[2][1] - &m[1][1] * &m[2][0]);
        if det != r(0) {
            matrix = candidate;
        }
    }
    let translation: [R; 3] = std::array::from_fn(|_| q(bytes.small(9, 4), 3));
    let transform = |p: &[R; 3]| -> [R; 3] {
        std::array::from_fn(|i| &translation[i] + (0..3).map(|j| &matrix[i][j] * &p[j]).sum::<R>())
    };
    let common = pow2(i32::from(bytes.next() as i8) * 8);
    let controls = case
        .controls
        .iter()
        .map(|c| {
            let w = &c[3];
            let [x, y, z] = std::array::from_fn(|i| {
                &translation[i] * w + (0..3).map(|j| &matrix[i][j] * &c[j]).sum::<R>()
            });
            [x * &common, y * &common, z * &common, w * &common]
        })
        .collect();
    let (a, b) = (transform(&case.a), transform(&case.b));
    let knots = case.knots.iter().map(global).collect();
    let basis = if case.periodic {
        ExactKnotVector::new_periodic(case.degree, knots, case.mults.clone())
    } else {
        ExactKnotVector::new(case.degree, knots, case.mults.clone())
    }
    .unwrap();
    let mut curve = ExactBSplineCurve3::from_homogeneous(basis, controls).unwrap();
    // Exact edits preserve the parameter function and hence the answer.
    let edit = bytes.next();
    let [d0, d1] = curve.domain().clone();
    if edit & 1 != 0 {
        let cut = &d0 + (&d1 - &d0) * q(1 + i64::from(bytes.next()), 257);
        curve = curve.insert_knot(&cut, case.degree).unwrap();
    }
    if edit & 2 != 0 && case.degree < 25 {
        let target = (case.degree + 1 + usize::from(bytes.next() % 25)).min(25);
        curve = curve.elevated(target).unwrap();
    }
    let (first, last) = (global(&case.lo), global(&case.hi));
    let run = |a: &[R; 3], b: &[R; 3]| {
        if line {
            exact_spline_line_in(&curve, a, b, &first, &last)
        } else {
            exact_spline_segment_in(&curve, a, b, &first, &last)
        }
        .unwrap()
    };
    let result = run(&a, &b);
    assert_eq!(result.points().len(), case.points.len());
    assert_eq!(result.overlaps().len(), case.intervals.len());
    assert_eq!(
        result.is_disjoint(),
        case.points.is_empty() && case.intervals.is_empty()
    );
    let verify = |actual: &SplineLinearPoint, hit: &Hit| {
        identity("parameter", &hit.t, &offset, &width, |v| {
            actual.compare_parameter(v).unwrap()
        });
        enclosure(actual.parameter_bounds(), &hit.t.affine(&offset, &width));
        let coordinates: [Value; 3] = std::array::from_fn(|i| {
            let k: R = (0..3).map(|j| &matrix[i][j] * &hit.direction[j]).sum();
            identity("coordinate", &hit.scale, &translation[i], &k, |v| {
                actual.compare_coordinate(i, v).unwrap()
            });
            hit.scale.affine(&translation[i], &k)
        });
        if coordinates.iter().all(representable) {
            let actual = actual.coordinate_bounds().unwrap();
            for (bound, expected) in actual.into_iter().zip(&coordinates) {
                enclosure(Ok(bound), expected);
            }
        } else {
            assert!(matches!(
                actual.coordinate_bounds(),
                Err(Error::Unrepresentable(_))
            ));
        }
        identity("linear parameter", &hit.linear, &r(0), &r(1), |v| {
            actual.compare_linear_parameter(v).unwrap()
        });
        enclosure(actual.linear_parameter_bounds(), &hit.linear);
    };
    for (actual, hit) in result.points().iter().zip(&case.points) {
        verify(actual, hit);
    }
    for (actual, hits) in result.overlaps().iter().zip(&case.intervals) {
        for (endpoint, hit) in actual.endpoints().iter().zip(hits) {
            verify(endpoint, hit);
        }
        assert_eq!(
            actual.endpoints()[0].parameter_cmp(&actual.endpoints()[1]),
            Less
        );
    }
    // Reversing the query keeps every curve parameter and maps s to 1-s.
    let reversed = run(&b, &a);
    assert_eq!(reversed.points().len(), case.points.len());
    assert_eq!(reversed.overlaps().len(), case.intervals.len());
    let flipped = |actual: &SplineLinearPoint, expected: &SplineLinearPoint, hit: &Hit| {
        assert_eq!(actual.parameter_cmp(expected), Equal);
        let s = if a == b {
            hit.linear.clone()
        } else {
            hit.linear.affine(&r(1), &-r(1))
        };
        identity("reversed linear parameter", &s, &r(0), &r(1), |v| {
            actual.compare_linear_parameter(v).unwrap()
        });
    };
    for ((x, y), hit) in reversed
        .points()
        .iter()
        .zip(result.points())
        .zip(&case.points)
    {
        flipped(x, y, hit);
    }
    for ((x, y), hits) in reversed
        .overlaps()
        .iter()
        .zip(result.overlaps())
        .zip(&case.intervals)
    {
        for ((p, e), hit) in x.endpoints().iter().zip(y.endpoints()).zip(hits) {
            flipped(p, e, hit);
        }
    }
    // Typed rejections and atomic work limits.
    let bad = R::new_raw(1.into(), 0.into());
    let malformed = [bad.clone(), r(0), r(0)];
    assert!(exact_spline_segment_in(&curve, &malformed, &b, &first, &last).is_err());
    assert!(exact_spline_line_in(&curve, &a, &a, &first, &last).is_err());
    if first < last {
        assert!(matches!(
            exact_spline_segment_in(&curve, &a, &b, &last, &first),
            Err(Error::OutOfDomain(_))
        ));
    }
    assert!(matches!(
        exact_spline_segment_in_with_options(
            &curve,
            &a,
            &b,
            &first,
            &last,
            SplineLinearOptions {
                max_spans: 0,
                ..Default::default()
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
    if let Some(p) = result.points().first() {
        assert!(p.compare_parameter(&bad).is_err());
        assert!(p.compare_linear_parameter(&bad).is_err());
        assert!(matches!(
            p.compare_coordinate(3, &r(0)),
            Err(Error::OutOfDomain(_))
        ));
    }
}
