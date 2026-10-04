//! S9d.4b.2a: a whole torus against prisms with arcs, spheres and cones
//! (REVIEW_NOTES.md, S9d.4b.2 split). On the torus's angles `(u, v)` the
//! other quadric's function `G(u, v)` is a trigonometric polynomial of degree
//! two in each: along the tube's circle at a rational `u` its roots are a
//! quartic's in the half-angle tangent of `v`, so the meeting's point there
//! is algebraic of degree four (`Q(alpha)`), and symmetrically at a rational
//! `v`. The meeting is cut into graphs over `u` where its slope in the angles
//! is at most one (away from its turning points in `u`, where `G_v = 0`) and
//! over `v` elsewhere (away from those in `v`), switched at rational
//! parameters where the slope is one: D13's procedural cell,
//! `Curve3::Toric`, each piece a graph over
//! its parameter's range whose other angle is the one root of `G` in a
//! rational window. Each piece is verified exactly: no root of `G` on the
//! window's ends over the range (Sturm counts), one root inside at one
//! parameter (isolation), and no double root inside the rectangle (certified
//! boxes, each clear of `G = 0` or of its derivative in the other angle).
//! The meeting's components are found by tracing from the roots at a
//! rational `u` in every interval between the critical values of `u` (the
//! roots of the discriminant of `G` in `v`); every root there must lie on a
//! verified piece, so no component is missed. At every critical value the
//! meeting must be regular (a turning point: `G_u != 0` wherever `G = G_v =
//! 0`, certified boxes): a singular point is a tangency of the surfaces,
//! `Degenerate`. Coaxial pairs (`G` independent of `u`) meet in circles:
//! rings over `u` at the roots in `v`. A cap's or rim's circle meets the
//! torus where the torus's quartic along it vanishes: a polynomial of degree
//! eight in its half-angle tangent.
//!
//! S9d.4b.2b: the other surface another whole torus. Its function `(|l|^2 +
//! R2^2 - r2^2)^2 - 4 R2^2 (l_u^2 + l_v^2)` in its own local coordinates `l`
//! is of degree two in each angle where both frames are exactly orthonormal
//! (a round circle keeps `|l|^2` affine in its angle's cosine and sine), as
//! a quadric's, and the meeting is traced the same way; in stored frames
//! not exactly orthonormal (rounded axes) it is of degree four in each, a
//! point at a rational parameter algebraic of degree eight, whose
//! discriminant is out of reach: its critical values are enclosed instead by
//! a certified subdivision of the angles (boxes clear of `G` or `G_v`, the
//! rest small and clear of `G_u`: turning points, a box clear of none at
//! `1e-10` a tangency), lines between them seed the traces, and every
//! turning point's box must lie in a verified piece over `v` (so no
//! component is missed).
use super::graph::{between_ccw, same_dir};
use super::meet::{CylPair, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::procedural::Other;
use super::torus::Ring;
use super::torus_segment::Span;
use super::turned::{
    changes, padd, pmul, pscale, roots, square_sum, sturm, trim, Chart, Form, Poly,
};
use crate::certified::{Fast, Real};
use crate::polynomial::real::{isolate, AlgebraicRoot, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::sync::{Arc, Mutex};

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

/// A singular point of the meeting: the surfaces tangent there.
fn tangent() -> Error {
    Error::Degenerate("a torus tangent to the other input's surface (S9d.4b.2)")
}

// ------------------------------------------------------------ forms

/// Exponents of `cos u`, `sin u`, `cos v`, `sin v`.
type Exp = [u32; 4];

/// A trigonometric polynomial on the torus's angles, `sum k cu^a su^b cv^c
/// sv^d`, reduced by `su^2 = 1 - cu^2` and `sv^2 = 1 - cv^2` (so `b, d <=
/// 1`): zero exactly when every coefficient is.
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct Bi {
    terms: BTreeMap<Exp, R>,
}

impl Bi {
    fn add_term(&mut self, e: Exp, k: &R) {
        if *k == zero() {
            return;
        }
        let t = self.terms.entry(e).or_insert_with(zero);
        *t += k;
        if *t == zero() {
            self.terms.remove(&e);
        }
    }

    fn add(&self, o: &Self) -> Self {
        let mut out = self.clone();
        for (e, k) in &o.terms {
            out.add_term(*e, k);
        }
        out
    }

    fn scale(&self, k: &R) -> Self {
        let mut out = Self::default();
        for (e, x) in &self.terms {
            out.add_term(*e, &(x * k));
        }
        out
    }

    fn sub(&self, o: &Self) -> Self {
        self.add(&o.scale(&int(-1)))
    }

    fn mul(&self, o: &Self) -> Self {
        let mut out = Self::default();
        for (a, x) in &self.terms {
            for (b, y) in &o.terms {
                let e = [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
                out.add_term(e, &(x * y));
            }
        }
        out.reduce()
    }

    /// `su^2 = 1 - cu^2`, `sv^2 = 1 - cv^2` until every `b, d <= 1`.
    fn reduce(self) -> Self {
        let mut out = Self::default();
        let mut stack: Vec<(Exp, R)> = self.terms.into_iter().collect();
        while let Some((e, k)) = stack.pop() {
            if e[1] >= 2 {
                stack.push(([e[0], e[1] - 2, e[2], e[3]], k.clone()));
                stack.push(([e[0] + 2, e[1] - 2, e[2], e[3]], -k));
            } else if e[3] >= 2 {
                stack.push(([e[0], e[1], e[2], e[3] - 2], k.clone()));
                stack.push(([e[0], e[1], e[2] + 2, e[3] - 2], -k));
            } else {
                out.add_term(e, &k);
            }
        }
        out
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Whether it depends on `u` (false: coaxial).
    fn depends_on_u(&self) -> bool {
        self.terms.keys().any(|e| e[0] > 0 || e[1] > 0)
    }

    /// The derivative in `u` (`of_v` false) or in `v`.
    fn d(&self, of_v: bool) -> Self {
        let (i, j) = if of_v { (2, 3) } else { (0, 1) };
        let mut out = Self::default();
        for (e, k) in &self.terms {
            // d cos^a sin^b = -a cos^(a-1) sin^(b+1) + b cos^(a+1) sin^(b-1).
            if e[i] > 0 {
                let mut f = *e;
                f[i] -= 1;
                f[j] += 1;
                out.add_term(f, &(k * int(-i64::from(e[i]))));
            }
            if e[j] > 0 {
                let mut f = *e;
                f[i] += 1;
                f[j] -= 1;
                out.add_term(f, &(k * int(i64::from(e[j]))));
            }
        }
        out.reduce()
    }

    /// Its degree in `u` (`of_v` false) or in `v`.
    fn degree(&self, of_v: bool) -> u32 {
        let (i, j) = if of_v { (2, 3) } else { (0, 1) };
        self.terms.keys().map(|e| e[i] + e[j]).max().unwrap_or(0)
    }

    /// The form in the other angle at a rational `(cos, sin)` of `u`
    /// (`of_v` false) or of `v`: of degree two (four for two tori in frames
    /// not exactly orthonormal).
    fn at(&self, of_v: bool, cs: &[R; 2]) -> Form {
        let pow = |x: &R, n: u32| (0..n).fold(int(1), |acc, _| acc * x);
        let (i, j, k, l) = if of_v { (2, 3, 0, 1) } else { (0, 1, 2, 3) };
        let mut terms: BTreeMap<(u32, u32), R> = BTreeMap::new();
        for (e, x) in &self.terms {
            let c = x * pow(&cs[0], e[i]) * pow(&cs[1], e[j]);
            *terms.entry((e[k], e[l])).or_insert_with(zero) += c;
        }
        Form::from_terms(terms, self.degree(!of_v).max(2))
    }

    /// The form in the other angle at a `(cos, sin)` of `u` (`of_v` false)
    /// or of `v` in one quadratic field `Q(sqrt(d))` (S9d.4c: a rim's fixed
    /// angle): `A + sqrt(d) B` with `A`, `B` rational forms, and `d` (zero
    /// where the direction is rational).
    pub(super) fn at_dir(&self, of_v: bool, cs: &[Qd; 2]) -> (Form, Form, R) {
        let d = cs
            .iter()
            .find_map(|c| c.field().cloned())
            .unwrap_or_else(zero);
        let pow = |x: &Qd, n: u32| (0..n).fold(Qd::rat(int(1)), |acc, _| acc.mul(x));
        let (i, j, k, l) = if of_v { (2, 3, 0, 1) } else { (0, 1, 2, 3) };
        let mut ta: BTreeMap<(u32, u32), R> = BTreeMap::new();
        let mut tb: BTreeMap<(u32, u32), R> = BTreeMap::new();
        for (e, x) in &self.terms {
            let c = pow(&cs[0], e[i]).mul(&pow(&cs[1], e[j]));
            let part = |v: &K| match v {
                K::Rat(r) => r.clone(),
                _ => unreachable!("a direction of Q(sqrt(d))"),
            };
            let (a, b) = if c.field().is_some() {
                (part(&c.a), part(&c.b))
            } else {
                (part(&c.a), zero())
            };
            *ta.entry((e[k], e[l])).or_insert_with(zero) += x * a;
            *tb.entry((e[k], e[l])).or_insert_with(zero) += x * b;
        }
        let deg = self.degree(!of_v).max(2);
        (Form::from_terms(ta, deg), Form::from_terms(tb, deg), d)
    }

    fn to_f64(&self) -> F64Bi {
        F64Bi {
            terms: self
                .terms
                .iter()
                .map(|(e, k)| (*e, rational_f64(k)))
                .collect(),
        }
    }

    fn to_fast(&self) -> FastBi {
        FastBi {
            terms: self
                .terms
                .iter()
                .map(|(e, k)| (*e, Fast::from_r(k)))
                .collect(),
        }
    }
}

/// A binary64 view of a form.
#[derive(Debug, Clone)]
struct F64Bi {
    terms: Vec<(Exp, f64)>,
}

impl F64Bi {
    /// Its value at `(cu, su, cv, sv)`.
    fn eval(&self, c: &[f64; 4]) -> f64 {
        let pow = |x: f64, n: u32| match n {
            0 => 1.0,
            1 => x,
            2 => x * x,
            _ => x.powi(n as i32),
        };
        self.terms
            .iter()
            .map(|(e, k)| k * pow(c[0], e[0]) * pow(c[1], e[1]) * pow(c[2], e[2]) * pow(c[3], e[3]))
            .sum()
    }
}

/// A certified binary64 view of a form (its coefficients enclosed).
#[derive(Debug, Clone)]
struct FastBi {
    terms: Vec<(Exp, Fast)>,
}

impl FastBi {
    /// An enclosure over enclosed `(cu, su, cv, sv)` (natural extension).
    fn eval(&self, c: &[Fast; 4]) -> Fast {
        let pow = |x: &Fast, n: u32| match n {
            0 => Fast::exact_f64(1.0),
            1 => *x,
            2 => x.square(),
            3 => x.square().mul(x),
            4 => x.square().square(),
            _ => (0..n).fold(Fast::exact_f64(1.0), |acc, _| acc.mul(x)),
        };
        self.terms.iter().fold(Fast::exact_f64(0.0), |acc, (e, k)| {
            acc.add(
                &k.mul(&pow(&c[0], e[0]))
                    .mul(&pow(&c[1], e[1]))
                    .mul(&pow(&c[2], e[2]))
                    .mul(&pow(&c[3], e[3])),
            )
        })
    }
}

/// The cosines and sines of a box's angles, enclosed.
fn trig_box(u: &Fast, v: &Fast) -> [Fast; 4] {
    let (cu, su) = Fast::cos_sin(u);
    let (cv, sv) = Fast::cos_sin(v);
    [cu, su, cv, sv]
}

fn fast_span(lo: f64, hi: f64) -> Fast {
    Fast::exact_f64(lo).union(&Fast::exact_f64(hi))
}

/// A function and its two derivatives, certified: its range over a box by
/// the natural extension and the mean-value form, intersected.
#[derive(Debug, Clone)]
struct Graded {
    f: FastBi,
    fu: FastBi,
    fv: FastBi,
}

impl Graded {
    fn new(f: &Bi) -> Self {
        Self {
            f: f.to_fast(),
            fu: f.d(false).to_fast(),
            fv: f.d(true).to_fast(),
        }
    }

    /// Whether the function is certainly nonzero over the box `[u] x [v]`.
    fn clear(&self, u: [f64; 2], v: [f64; 2]) -> bool {
        let whole = trig_box(&fast_span(u[0], u[1]), &fast_span(v[0], v[1]));
        let natural = self.f.eval(&whole);
        if natural.sign().is_some_and(|s| s != Ordering::Equal) {
            return true;
        }
        let (mu, mv) = (0.5 * u[0] + 0.5 * u[1], 0.5 * v[0] + 0.5 * v[1]);
        let centre = self
            .f
            .eval(&trig_box(&Fast::exact_f64(mu), &Fast::exact_f64(mv)));
        let du = fast_span(u[0], u[1]).sub(&Fast::exact_f64(mu));
        let dv = fast_span(v[0], v[1]).sub(&Fast::exact_f64(mv));
        let mean = centre
            .add(&self.fu.eval(&whole).mul(&du))
            .add(&self.fv.eval(&whole).mul(&dv));
        mean.sign().is_some_and(|s| s != Ordering::Equal)
    }
}

// ------------------------------------------------------------ the meeting

/// The other input's surface: a quadric (S9d.4b.2a) or a whole torus
/// (S9d.4b.2b, its model's frame and radii).
#[derive(Debug, Clone)]
pub(super) enum Far {
    Quadric(Box<Other>),
    Torus { f: Box<Affine>, big: R, small: R },
}

impl Far {
    fn ring(big: &R, small: &R) -> Ring {
        Ring {
            big: big.clone(),
            small: small.clone(),
            e: [int(1), zero()],
            v0: [int(1), zero()],
            span: Span::Whole,
            reversed: false,
            rims: [None, None],
        }
    }

    /// Its function's exact value at a point (zero on it).
    fn value(&self, x: &QV) -> Qd {
        match self {
            Far::Quadric(o) => o.value(x),
            Far::Torus { f, big, small } => Self::ring(big, small).value(&f.local_q(x)),
        }
    }

    /// A positive multiple of its function's gradient at a point.
    fn gradient(&self, x: &QV) -> QV {
        match self {
            Far::Quadric(o) => o.gradient(x),
            Far::Torus { f, big, small } => {
                let g = Self::ring(big, small).gradient(&f.local_q(x));
                qadd(
                    &qadd(&qscale(f.row(0), &g[0]), &qscale(f.row(1), &g[1])),
                    &qscale(f.row(2), &g[2]),
                )
            }
        }
    }
}

/// The other surface's function on a torus's angles, `G(u, v)`: its
/// quadric's (S9d.4b.2a) or its torus's (S9d.4b.2b) at the torus's point
/// `o + (R + r cv)(cu x + su y) + r sv n`, reduced.
pub(super) fn meeting_form(f: &Affine, big: &R, small: &R, other: &Far) -> Bi {
    // The torus's point: o + (R + r cv)(cu x + su y) + r sv n.
    let lin = |g: &V, e: &R| -> Bi {
        let mut b = Bi::default();
        b.add_term([0, 0, 0, 0], &(dot(g, &f.o) - e));
        b.add_term([1, 0, 0, 0], &(big * dot(g, &f.x)));
        b.add_term([0, 1, 0, 0], &(big * dot(g, &f.y)));
        b.add_term([1, 0, 1, 0], &(small * dot(g, &f.x)));
        b.add_term([0, 1, 1, 0], &(small * dot(g, &f.y)));
        b.add_term([0, 0, 0, 1], &(small * dot(g, &f.n)));
        b
    };
    match other {
        Far::Quadric(other) => {
            let mut g = Bi::default();
            for (gi, ei) in other.g.iter().zip(&other.e) {
                let l = lin(gi, ei);
                g = g.add(&l.mul(&l));
            }
            // The radius term r + t (h . p - e_h).
            let mut rad = lin(&other.h, &other.eh).scale(&other.t);
            rad.add_term([0, 0, 0, 0], &other.r);
            g.sub(&rad.mul(&rad)).reduce()
        }
        Far::Torus {
            f: f2,
            big: b2,
            small: s2,
        } => {
            // Its local coordinates `l_k = row_k . (p - o2)`: `S^2 - 4
            // R2^2 P`, `S = |l|^2 + R2^2 - r2^2`, `P = l_u^2 + l_v^2`.
            let l: Vec<Bi> = (0..3)
                .map(|k| lin(f2.row(k), &dot(f2.row(k), &f2.o)))
                .collect();
            let p = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
            let mut sum = p.add(&l[2].mul(&l[2]));
            sum.add_term([0, 0, 0, 0], &(b2 * b2 - s2 * s2));
            sum.mul(&sum).sub(&p.scale(&(int(4) * b2 * b2))).reduce()
        }
    }
}

/// A torus and a quadric or another torus (the other input's face): `G`
/// on the torus's angles, exact and in binary64, and its certified views.
#[derive(Debug)]
pub(super) struct Meeting {
    /// The torus's frame and radii (its model's).
    pub(super) f: Affine,
    pub(super) big: R,
    pub(super) small: R,
    pub(super) other: Far,
    g: Bi,
    num: [F64Bi; 3],
    /// `G` and its derivative in `u` and in `v`, certified with their own
    /// derivatives.
    cert: [Graded; 3],
    /// The binary64 frame, for samples.
    frame: [[f64; 3]; 4],
    /// The sum of `G`'s coefficients' sizes (its binary64 values' scale).
    scale: f64,
    /// Tangents asked for, by point, with their runs' signs in `u` and `v`.
    tangents: Mutex<Vec<(QV, QV, [Ordering; 2])>>,
}

impl Meeting {
    pub(super) fn new(f: &Affine, big: &R, small: &R, other: &Far) -> Self {
        let g = meeting_form(f, big, small, other);
        let (gu, gv) = (g.d(false), g.d(true));
        let fl = |v: &V| v.clone().map(|x| rational_f64(&x));
        let num = [g.to_f64(), gu.to_f64(), gv.to_f64()];
        let scale = num[0]
            .terms
            .iter()
            .map(|(_, k)| k.abs())
            .sum::<f64>()
            .max(1.0);
        Self {
            f: f.clone(),
            big: big.clone(),
            small: small.clone(),
            other: other.clone(),
            num,
            cert: [Graded::new(&g), Graded::new(&gu), Graded::new(&gv)],
            g,
            frame: [fl(&f.o), fl(&f.x), fl(&f.y), fl(&f.n)],
            scale,
            tangents: Mutex::new(Vec::new()),
        }
    }

    fn ring(&self) -> Ring {
        Far::ring(&self.big, &self.small)
    }

    /// The torus's and the other surface's gradients crossed at a point,
    /// and the signs of its runs in `u` and in `v`.
    fn tangent_runs(&self, x: &QV) -> (QV, [Ordering; 2]) {
        let f = &self.f;
        let l = f.local_q(x);
        let ring = self.ring();
        let g = ring.gradient(&l);
        let row = |k: usize| f.row(k).clone();
        let gt = qadd(
            &qadd(&qscale(&row(0), &g[0]), &qscale(&row(1), &g[1])),
            &qscale(&row(2), &g[2]),
        );
        let t = qcross(&gt, &self.other.gradient(x));
        let lt = f.local_dir_q(&t);
        let run_u = l[0].mul(&lt[1]).sub(&l[1].mul(&lt[0])).sign();
        let rho = ring.rho(&l);
        let drho = l[0].mul(&lt[0]).add(&l[1].mul(&lt[1]));
        let run_v = rho
            .add_r(&-self.big.clone())
            .mul(&rho)
            .mul(&lt[2])
            .sub(&l[2].mul(&drho))
            .sign();
        (t, [run_u, run_v])
    }

    /// `G`, `G_u` and `G_v` in binary64.
    fn eval(&self, u: f64, v: f64) -> [f64; 3] {
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let c = [cu, su, cv, sv];
        [
            self.num[0].eval(&c),
            self.num[1].eval(&c),
            self.num[2].eval(&c),
        ]
    }

    /// `G` at a parameter `t` (`u`, or `v` when `over_v`) and the other
    /// angle `s`, in binary64.
    fn value(&self, over_v: bool, t: f64, s: f64) -> f64 {
        let (u, v) = if over_v { (s, t) } else { (t, s) };
        self.eval(u, v)[0]
    }

    /// The torus's binary64 point at its angles.
    fn point_f64(&self, u: f64, v: f64) -> [f64; 3] {
        let [o, x, y, n] = self.frame;
        let (big, small) = (rational_f64(&self.big), rational_f64(&self.small));
        let rho = big + small * v.cos();
        let l = [rho * u.cos(), rho * u.sin(), small * v.sin()];
        [0, 1, 2].map(|j| o[j] + l[0] * x[j] + l[1] * y[j] + l[2] * n[j])
    }

    /// The torus's exact point at `(cos, sin)` of `u` and of `v`.
    fn point_q(&self, cu: &[Qd; 2], cv: &[Qd; 2]) -> QV {
        let rho = cv[0].scale(&self.small).add_r(&self.big);
        let l = [rho.mul(&cu[0]), rho.mul(&cu[1]), cv[1].scale(&self.small)];
        qadd(
            &qv(&self.f.o),
            &qadd(
                &qadd(&qscale(&self.f.x, &l[0]), &qscale(&self.f.y, &l[1])),
                &qscale(&self.f.n, &l[2]),
            ),
        )
    }

    /// A point's direction of `u` (`of_v` false) or of `v`, exact (on the
    /// torus): a positive multiple of its `(cos, sin)`, `(l_u, l_v)` or
    /// `(rho - R, l_w)`.
    fn dir_raw(&self, of_v: bool, x: &QV) -> [Qd; 2] {
        let l = self.f.local_q(x);
        if of_v {
            self.ring().v_dir(&l)
        } else {
            Ring::u_dir(&l)
        }
    }

    /// The roots in the other angle at a rational parameter, exactly: its
    /// form's polynomial in a chart whose antipode is no root, and the
    /// roots (a repeated one refused).
    fn roots_at(&self, over_v: bool, p: &[R; 2]) -> Result<(Chart, Poly, Vec<AlgebraicRoot>)> {
        let form = self.g.at(over_v, p);
        let chart = clear_chart(&form)?;
        let poly = trim(form.poly(&chart));
        if poly.is_empty() {
            return Err(tangent());
        }
        let rs = roots(&poly)?;
        Ok((chart, poly, rs))
    }

    /// The exact point at a rational parameter `p` whose other angle is the
    /// root `root` of `poly` in `chart`, with its angles' directions.
    fn exact_point(
        &self,
        over_v: bool,
        p: &[R; 2],
        chart: &Chart,
        poly: &Poly,
        root: AlgebraicRoot,
    ) -> Result<MeetPoint> {
        let g = Arc::new(Gen::new(poly.clone(), root));
        let t = K::generator(&g);
        let den = t.mul(&t).add(&K::Rat(int(1)));
        // `(1 - t^2, 2 t)` turned by the chart's base: a positive multiple
        // of the other angle's `(cos, sin)`.
        let (c, s) = (K::Rat(int(1)).sub(&t.mul(&t)), t.scale(&int(2)));
        let raw = [
            c.scale(&chart.c0).sub(&s.scale(&chart.s0)),
            c.scale(&chart.s0).add(&s.scale(&chart.c0)),
        ];
        let inv = den
            .recip()
            .ok_or(limit("a torus meeting's point at infinity"))?;
        let other = raw.clone().map(|x| Qd::of(x.mul(&inv)));
        let param = [Qd::rat(p[0].clone()), Qd::rat(p[1].clone())];
        let raw = raw.map(Qd::of);
        Ok(if over_v {
            MeetPoint {
                x: self.point_q(&other, &param),
                u: raw,
                v: param,
            }
        } else {
            MeetPoint {
                x: self.point_q(&param, &other),
                u: param,
                v: raw,
            }
        })
    }
}

/// A point of the meeting with the directions of its angles `u` and `v`
/// (positive multiples of their `(cos, sin)`).
#[derive(Debug, Clone)]
struct MeetPoint {
    x: QV,
    u: [Qd; 2],
    v: [Qd; 2],
}

/// A chart of the circle whose antipode the form does not vanish at.
fn clear_chart(form: &Form) -> Result<Chart> {
    for (c0, s0) in [(1, 0), (0, 1), (-1, 0), (0, -1)] {
        if form.value(&[int(-c0), int(-s0)]) != zero() {
            return Ok(Chart {
                c0: int(c0),
                s0: int(s0),
            });
        }
    }
    for k in 1..40i64 {
        let d = circle_point(&[zero(), zero()], &int(1), &(int(k) / int(7)));
        if form.value(&[-d[0].clone(), -d[1].clone()]) != zero() {
            return Ok(Chart {
                c0: d[0].clone(),
                s0: d[1].clone(),
            });
        }
    }
    Err(tangent())
}

/// A rational unit direction near a binary64 angle.
fn dir_near(angle: f64) -> [R; 2] {
    let half = angle / 2.0;
    let one = int(1);
    if half.cos().abs() >= 0.5 {
        let s = q(half.tan());
        let den = &one + &s * &s;
        [(&one - &s * &s) / &den, int(2) * &s / &den]
    } else {
        let c = q(half.cos() / half.sin());
        let den = &c * &c + &one;
        [(&c * &c - &one) / &den, int(2) * &c / &den]
    }
}

fn angle_r(d: &[R; 2]) -> f64 {
    rational_f64(&d[1]).atan2(rational_f64(&d[0]))
}

fn angle_q(d: &[Qd; 2]) -> f64 {
    d[1].to_f64().atan2(d[0].to_f64())
}

fn qdir(d: &[R; 2]) -> [Qd; 2] {
    [Qd::rat(d[0].clone()), Qd::rat(d[1].clone())]
}

/// `x` lifted by whole turns nearest `target`.
fn near(x: f64, target: f64) -> f64 {
    x + TAU * ((target - x) / TAU).round()
}

/// The binary64 angle of a root of a chart's polynomial.
fn root_angle(chart: &Chart, r: &AlgebraicRoot) -> f64 {
    let mut r = r.clone();
    r.refine_for_signs(72);
    let (a, b) = r.isolator();
    angle_r(&chart.at(&((a + b) / int(2))))
}

/// A polynomial's exact value at a number of any field.
fn peval(p: &Poly, x: &Qd) -> Qd {
    let mut acc = Qd::rat(zero());
    for c in p.iter().rev() {
        acc = acc.mul(x).add_r(c);
    }
    acc
}

/// A root in `[a, b]` of a binary64 function changing sign there (bisection
/// with secant steps well inside).
fn bracket_root(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let (mut a, mut b) = (a, b);
    let (mut fa, mut fb) = (f(a), f(b));
    if fa.signum() == fb.signum() && fa != 0.0 && fb != 0.0 {
        let n = 64;
        let mut prev = (a, fa);
        for k in 1..=n {
            let x = a + (b - a) * k as f64 / n as f64;
            let fx = f(x);
            if fx.signum() != prev.1.signum() {
                (a, fa, b, fb) = (prev.0, prev.1, x, fx);
                break;
            }
            prev = (x, fx);
        }
    }
    if fa == 0.0 {
        return a;
    }
    if fb == 0.0 {
        return b;
    }
    for _ in 0..200 {
        if b - a <= 4.0 * f64::EPSILON * a.abs().max(b.abs()).max(1.0) {
            break;
        }
        let mid = 0.5 * a + 0.5 * b;
        let secant = a - fa * (b - a) / (fb - fa);
        let x = if secant > a && secant < b && (secant - mid).abs() < 0.25 * (b - a) {
            secant
        } else {
            mid
        };
        let fx = f(x);
        if fx == 0.0 {
            return x;
        }
        if fx.signum() == fa.signum() {
            (a, fa) = (x, fx);
        } else {
            (b, fb) = (x, fx);
        }
    }
    0.5 * a + 0.5 * b
}

// ------------------------------------------------------------ pieces

/// A piece of a torus's meeting with a quadric: a graph over `u` (or over
/// `v` when `over_v`) whose other angle is the one root of `G` strictly
/// within the counter-clockwise `window`, over the counter-clockwise `range`
/// of its parameter's direction (a ring: none).
#[derive(Debug, Clone)]
pub(super) struct ToricCrv {
    /// The torus's operand.
    pub(super) carrier: usize,
    pub(super) m: Arc<Meeting>,
    pub(super) over_v: bool,
    pub(super) window: [P2; 2],
    /// The window's middle direction: its chart's base.
    mid: P2,
    pub(super) range: Option<[[Qd; 2]; 2]>,
    /// The other quadric coaxial with the torus: a circle about its axis.
    pub(super) coaxial: bool,
    /// The window's and the range's binary64 angles, ascending.
    window_f: [f64; 2],
    range_f: Option<[f64; 2]>,
}

impl PartialEq for ToricCrv {
    fn eq(&self, o: &Self) -> bool {
        self.carrier == o.carrier
            && Arc::ptr_eq(&self.m, &o.m)
            && self.over_v == o.over_v
            && self.window == o.window
            && self.range == o.range
    }
}

impl ToricCrv {
    fn new(
        carrier: usize,
        m: &Arc<Meeting>,
        over_v: bool,
        window: [f64; 2],
        range: Option<[[Qd; 2]; 2]>,
    ) -> Result<Self> {
        let w = [dir_near(window[0]), dir_near(window[1])];
        let mid = dir_near(0.5 * window[0] + 0.5 * window[1]);
        if !between_ccw(&qdir(&w[0]), &qdir(&mid), &qdir(&w[1])) {
            return Err(limit("a torus meeting's window within rounding"));
        }
        let ascending = |a: f64, b: f64| {
            let mut b = b;
            while b <= a {
                b += TAU;
            }
            [a, b]
        };
        let window_f = ascending(angle_r(&w[0]), angle_r(&w[1]));
        let range_f = range
            .as_ref()
            .map(|[lo, hi]| ascending(angle_q(lo), angle_q(hi)));
        Ok(Self {
            carrier,
            m: m.clone(),
            over_v,
            window: w,
            mid,
            range,
            coaxial: !m.g.depends_on_u(),
            window_f,
            range_f,
        })
    }

    fn chart(&self) -> Chart {
        Chart {
            c0: self.mid[0].clone(),
            s0: self.mid[1].clone(),
        }
    }

    /// The window's ends in its chart (rational).
    fn window_t(&self) -> Option<[R; 2]> {
        let chart = self.chart();
        let t0 = chart.t_of(&qdir(&self.window[0]))?;
        let t1 = chart.t_of(&qdir(&self.window[1]))?;
        Some([t0.rational()?.clone(), t1.rational()?.clone()])
    }

    /// The window's binary64 angles, ascending.
    pub(super) fn window_angles(&self) -> [f64; 2] {
        let a = angle_r(&self.window[0]);
        let mut b = angle_r(&self.window[1]);
        while b <= a {
            b += TAU;
        }
        [a, b]
    }

    /// The roots strictly within the window at a rational parameter, with
    /// their chart and polynomial.
    fn window_roots(&self, cs: &[R; 2]) -> Option<(Chart, Poly, Vec<AlgebraicRoot>)> {
        let chart = self.chart();
        let poly = trim(self.m.g.at(self.over_v, cs).poly(&chart));
        if poly.is_empty() {
            return None;
        }
        let ip = IntPolynomial::from_rationals(&poly);
        if ip.is_constant() {
            return Some((chart, poly, Vec::new()));
        }
        let [t0, t1] = self.window_t()?;
        let rs = isolate(
            &ip,
            t0.clone(),
            t1.clone(),
            &mut Budget::new(RootIsolationOptions::default()),
        )
        .ok()?
        .into_iter()
        .filter(|r| {
            r.compare_rational(&t0) == Ordering::Greater
                && r.compare_rational(&t1) == Ordering::Less
        })
        .collect();
        Some((chart, poly, rs))
    }

    /// The piece's point at a rational `(cos, sin)` of its parameter
    /// (`None` where the window holds no root, or more than one).
    pub(super) fn at(&self, cs: &[R; 2]) -> Option<QV> {
        Some(self.meet_at(cs)?.x)
    }

    /// The same with its angles' directions.
    fn meet_at(&self, cs: &[R; 2]) -> Option<MeetPoint> {
        let (chart, poly, mut rs) = self.window_roots(cs)?;
        if rs.len() != 1 {
            return None;
        }
        self.m
            .exact_point(self.over_v, cs, &chart, &poly, rs.pop()?)
            .ok()
    }

    /// Whether a point may lie on the piece: its binary64 view on both
    /// surfaces, its angles within the window and the range, each within a
    /// margin far above the view's error (a filter before the exact tests).
    fn near_f64(&self, x: &QV) -> bool {
        self.near_view(qv_f64(x))
    }

    /// `near_f64` at the point's binary64 view `p` (`qv_f64`).
    fn near_view(&self, p: [f64; 3]) -> bool {
        let [o, xa, ya, na] = self.m.frame;
        let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
        let l = super::graph::solve3(&xa, &ya, &na, &d);
        let (big, small) = (rational_f64(&self.m.big), rational_f64(&self.m.small));
        let rho = l[0].hypot(l[1]);
        if ((rho - big).hypot(l[2]) - small).abs() > 1e-7 * (1.0 + big) {
            return false;
        }
        let (u, v) = (l[1].atan2(l[0]), l[2].atan2(rho - big));
        let (t, s) = if self.over_v { (v, u) } else { (u, v) };
        let margin = 1e-7;
        let [w0, w1] = self.window_f;
        let s = near(s, 0.5 * w0 + 0.5 * w1);
        if s < w0 - margin || s > w1 + margin {
            return false;
        }
        if let Some([a, b]) = self.range_f {
            let t = near(t, 0.5 * a + 0.5 * b);
            if t < a - margin || t > b + margin {
                return false;
            }
        }
        let g = self.m.eval(u, v)[0];
        g.abs() <= 1e-7 * self.m.scale
    }

    /// A point's place: its parameter's direction (a positive multiple of
    /// its `(cos, sin)`).
    pub(super) fn place(&self, x: &QV) -> [Qd; 2] {
        self.m.dir_raw(self.over_v, x)
    }

    /// Whether a point of both surfaces lies on the piece: its other angle
    /// within the window, its parameter within the range (ends included).
    pub(super) fn holds(&self, x: &QV) -> bool {
        self.holds_at(x, qv_f64(x))
    }

    /// `holds` with the point's binary64 view (`qv_f64`) given: the
    /// arrangement tests each vertex against every piece of a meeting.
    pub(super) fn holds_at(&self, x: &QV, view: [f64; 3]) -> bool {
        if !self.near_view(view) {
            return false;
        }
        let (p, s) = (self.place(x), self.m.dir_raw(!self.over_v, x));
        self.holds_dirs(&p, &s)
    }

    /// Whether a meeting point lies on the piece, by its directions.
    fn holds_point(&self, x: &MeetPoint) -> bool {
        let (p, s) = if self.over_v {
            (&x.v, &x.u)
        } else {
            (&x.u, &x.v)
        };
        // Binary64 angles first, well clear of the window or the range.
        let margin = 1e-7;
        let [w0, w1] = self.window_f;
        let sa = near(angle_q(s), 0.5 * w0 + 0.5 * w1);
        if sa < w0 - margin || sa > w1 + margin {
            return false;
        }
        if let Some([a, b]) = self.range_f {
            let t = near(angle_q(p), 0.5 * a + 0.5 * b);
            if t < a - margin || t > b + margin {
                return false;
            }
        }
        self.holds_dirs(p, s)
    }

    /// Whether a parameter's and an other angle's directions lie within the
    /// range (ends included) and the window.
    fn holds_dirs(&self, p: &[Qd; 2], s: &[Qd; 2]) -> bool {
        if !between_ccw(&qdir(&self.window[0]), s, &qdir(&self.window[1])) {
            return false;
        }
        match &self.range {
            None => true,
            Some([lo, hi]) => same_dir(p, lo) || same_dir(p, hi) || between_ccw(lo, p, hi),
        }
    }

    /// Whether a point lies on the piece (both surfaces, the window, the
    /// range).
    pub(super) fn on(&self, x: &QV) -> bool {
        if !self.near_f64(x) {
            return false;
        }
        let l = self.m.f.local_q(x);
        self.m.ring().value(&l).sign() == Ordering::Equal
            && self.m.other.value(x).sign() == Ordering::Equal
            && self.holds(x)
    }

    /// The unit-free tangent at a point, running with the parameter: the
    /// torus's and the other surface's gradients crossed (kept by the
    /// meeting for each point asked, its pieces' and its edges' ends asked
    /// again and again).
    pub(super) fn tangent(&self, x: &QV) -> QV {
        let known = {
            let kept = self.m.tangents.lock().unwrap_or_else(|e| e.into_inner());
            kept.iter()
                .find(|(p, _, _)| p == x)
                .map(|(_, t, runs)| (t.clone(), runs[usize::from(self.over_v)]))
        };
        let (t, run) = match known {
            Some(k) => k,
            None => {
                let (t, runs) = self.m.tangent_runs(x);
                self.m
                    .tangents
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .push((x.clone(), t.clone(), runs));
                (t, runs[usize::from(self.over_v)])
            }
        };
        if run == Ordering::Less {
            t.map(|c| c.neg())
        } else {
            t
        }
    }

    /// The other angle at a binary64 parameter.
    pub(super) fn root_f64(&self, t: f64) -> f64 {
        let [a, b] = self.window_angles();
        bracket_root(|s| self.m.value(self.over_v, t, s), a, b)
    }

    /// Binary64 points from the parameter's angle `t0` turning `sweep`.
    pub(super) fn samples(&self, t0: f64, sweep: f64, n: usize) -> Vec<[f64; 3]> {
        let [a, b] = self.window_f;
        let mut prev: Option<f64> = None;
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                // Newton's steps from the last sample's root, kept inside
                // the window; the bracketing search otherwise.
                let from_prev = prev.and_then(|mut s| {
                    for _ in 0..8 {
                        let (u, v) = if self.over_v { (s, t) } else { (t, s) };
                        let [g, gu, gv] = self.m.eval(u, v);
                        let d = if self.over_v { gu } else { gv };
                        let step = g / d;
                        if !step.is_finite() {
                            return None;
                        }
                        s -= step;
                        if step.abs() <= 1e-15 * (1.0 + s.abs()) {
                            return (s > a && s < b).then_some(s);
                        }
                    }
                    None
                });
                let s = from_prev.unwrap_or_else(|| self.root_f64(t));
                prev = Some(s);
                if self.over_v {
                    self.m.point_f64(s, t)
                } else {
                    self.m.point_f64(t, s)
                }
            })
            .collect()
    }

    /// Whether the piece is verified: no root of `G` on the window's ends
    /// over the range (exact Sturm counts), one root inside at a parameter
    /// of the range (exact isolation), and no double root inside (certified
    /// boxes, each clear of `G` or of its derivative in the other angle).
    fn verified(&self) -> bool {
        // The window's ends.
        for w in &self.window {
            let form = self.m.g.at(!self.over_v, w);
            let clear = match &self.range {
                None => nowhere_zero(&form),
                Some([lo, hi]) => clear_over(&form, lo, hi),
            };
            if !clear {
                return false;
            }
        }
        // One root at a parameter inside.
        let probe = match &self.range {
            None => [int(1), zero()],
            Some([lo, hi]) => match rational_within(lo, hi) {
                Some(p) => p,
                None => return false,
            },
        };
        if self.at(&probe).is_none() {
            return false;
        }
        // No double root in the rectangle.
        let (p0, p1) = match &self.range {
            None => (0.0, TAU),
            Some([lo, hi]) => {
                let a = angle_q(lo);
                let mut b = angle_q(hi);
                while b <= a {
                    b += TAU;
                }
                (a - 1e-9, b + 1e-9)
            }
        };
        let [w0, w1] = self.window_angles();
        let g = &self.m.cert[0];
        let gs = &self.m.cert[if self.over_v { 1 } else { 2 }];
        let over_v = self.over_v;
        boxes_clear([p0, p1], [w0 - 1e-12, w1 + 1e-12], |p, s| {
            let (u, v) = if over_v { (s, p) } else { (p, s) };
            g.clear(u, v) || gs.clear(u, v)
        })
    }
}

/// Whether a form has no zero on the whole circle.
fn nowhere_zero(form: &Form) -> bool {
    if let Some(k) = form.as_constant() {
        return k != zero();
    }
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let poly = trim(form.poly(&chart));
    form.value(&[int(-1), zero()]) != zero()
        && !poly.is_empty()
        && matches!(super::turned::roots_repeated(&poly), Ok((_, r)) if r.is_empty())
}

/// Whether a form has no zero over the closed counter-clockwise range from
/// `lo` to `hi` (under a turn): exact Sturm counts in a chart based at its
/// middle.
fn clear_over(form: &Form, lo: &[Qd; 2], hi: &[Qd; 2]) -> bool {
    if let Some(k) = form.as_constant() {
        return k != zero();
    }
    let a = angle_q(lo);
    let mut b = angle_q(hi);
    while b <= a {
        b += TAU;
    }
    // Rational ends just outside algebraic ones (a zero between them is
    // refused as one inside: safe).
    let rational = |d: &[Qd; 2]| d[0].is_rational() && d[1].is_rational();
    if !rational(lo) || !rational(hi) {
        let (lo2, hi2) = (qdir(&dir_near(a - 1e-9)), qdir(&dir_near(b + 1e-9)));
        let inside = |x: &[Qd; 2]| between_ccw(&lo2, x, &hi2);
        if !inside(lo) || !inside(hi) || !between_ccw(&lo2, lo, hi) {
            return false;
        }
        return clear_over(form, &lo2, &hi2);
    }
    let mid = dir_near(0.5 * a + 0.5 * b);
    let chart = Chart {
        c0: mid[0].clone(),
        s0: mid[1].clone(),
    };
    let (Some(t0), Some(t1)) = (chart.t_of(lo), chart.t_of(hi)) else {
        return false;
    };
    if t0.cmp(&t1) != Ordering::Less {
        return false;
    }
    let poly = trim(form.poly(&chart));
    if poly.is_empty() {
        return false;
    }
    if peval(&poly, &t0).sign() == Ordering::Equal || peval(&poly, &t1).sign() == Ordering::Equal {
        return false;
    }
    if poly.len() == 1 {
        return true;
    }
    // The square-free part's chain (a repeated root counts once).
    let Ok((sf, _)) = super::turned::roots_repeated(&poly) else {
        return false;
    };
    // In integers at the rational ends (each member a positive multiple of
    // the rational chain's: the same changes).
    if let (Some(a), Some(b)) = (t0.rational(), t1.rational()) {
        let chain = super::turned::sturm_int(&sf);
        return super::turned::changes_int(&chain, a) == super::turned::changes_int(&chain, b);
    }
    let chain = sturm(&sf);
    changes(&chain, &t0) == changes(&chain, &t1)
}

/// A rational direction strictly within a counter-clockwise range.
fn rational_within(lo: &[Qd; 2], hi: &[Qd; 2]) -> Option<[R; 2]> {
    let a = angle_q(lo);
    let mut b = angle_q(hi);
    while b <= a {
        b += TAU;
    }
    for frac in [0.5, 0.25, 0.75, 0.375, 0.625] {
        let d = dir_near(a + frac * (b - a));
        let x = qdir(&d);
        if between_ccw(lo, &x, hi) {
            return Some(d);
        }
    }
    None
}

/// Whether every box of a subdivision of `[p] x [s]` is clear by `ok`
/// (halving the longer side, at most 60,000 boxes, none under 1e-10).
fn boxes_clear(p: [f64; 2], s: [f64; 2], ok: impl Fn([f64; 2], [f64; 2]) -> bool) -> bool {
    let mut stack = vec![(p, s)];
    let mut count = 0usize;
    while let Some((p, s)) = stack.pop() {
        count += 1;
        if count > 60_000 {
            return false;
        }
        if ok(p, s) {
            continue;
        }
        let (wp, ws) = (p[1] - p[0], s[1] - s[0]);
        if wp.max(ws) < 1e-10 {
            return false;
        }
        if wp >= ws {
            let m = 0.5 * p[0] + 0.5 * p[1];
            stack.push(([p[0], m], s));
            stack.push(([m, p[1]], s));
        } else {
            let m = 0.5 * s[0] + 0.5 * s[1];
            stack.push((p, [s[0], m]));
            stack.push((p, [m, s[1]]));
        }
    }
    true
}

// ------------------------------------------------------------ critical values

/// `G` times `(1 + t^2)^2 (1 + s^2)^2` in a chart of `u` (`t`) and one of
/// `v` (`s`): its coefficients of `s^j` as polynomials in `t`.
fn h_coefficients(g: &Bi, cu: &Chart, cv: &Chart) -> [Poly; 5] {
    let [cn, sn] = cu.numerators();
    let [cvn, svn] = cv.numerators();
    let w = vec![int(1), zero(), int(1)];
    let pw = |p: &Poly, n: u32| (0..n).fold(vec![int(1)], |acc, _| pmul(&acc, p));
    let mut out: [Poly; 5] = Default::default();
    for (e, k) in &g.terms {
        let rest_u = 2u32.saturating_sub(e[0] + e[1]);
        let rest_v = 2u32.saturating_sub(e[2] + e[3]);
        let pt = pscale(
            &pmul(&pmul(&pw(&cn, e[0]), &pw(&sn, e[1])), &pw(&w, rest_u)),
            k,
        );
        let ps = pmul(&pmul(&pw(&cvn, e[2]), &pw(&svn, e[3])), &pw(&w, rest_v));
        for (j, c) in ps.iter().enumerate() {
            out[j] = padd(&out[j], &pscale(&pt, c));
        }
    }
    out
}

/// The discriminant of the binary quartic `sum h_j s^j` (coefficients
/// polynomials in `t`).
/// In integers: the coefficients over their common denominator `D`, the
/// discriminant (of degree six in them) over `D^6`, each coefficient
/// reduced once: the same polynomial as in rationals reduced at every
/// product.
pub(super) fn quartic_discriminant(h: &[Poly; 5]) -> Poly {
    use num_bigint::BigInt;
    type IPoly = Vec<BigInt>;
    let den = h.iter().flatten().fold(BigInt::from(1), |m, x| {
        num_integer::Integer::lcm(&m, x.denom())
    });
    let scaled: Vec<IPoly> = h
        .iter()
        .map(|p| p.iter().map(|x| x.numer() * (&den / x.denom())).collect())
        .collect();
    let mul = |p: &IPoly, q: &IPoly| -> IPoly {
        if p.is_empty() || q.is_empty() {
            return Vec::new();
        }
        let mut out = vec![BigInt::from(0); p.len() + q.len() - 1];
        for (i, x) in p.iter().enumerate() {
            for (j, y) in q.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    };
    let (a, b, c, d, e) = (&scaled[4], &scaled[3], &scaled[2], &scaled[1], &scaled[0]);
    // Products shared by the terms.
    let (aa, bb, cc, dd, ee) = (mul(a, a), mul(b, b), mul(c, c), mul(d, d), mul(e, e));
    let (ae, bd) = (mul(a, e), mul(b, d));
    let terms: [(i64, Vec<&IPoly>); 16] = [
        (256, vec![&aa, a, e, &ee]),
        (-192, vec![&aa, &bd, &ee]),
        (-128, vec![&aa, &cc, &ee]),
        (144, vec![&aa, c, &dd, e]),
        (-27, vec![&aa, &dd, &dd]),
        (144, vec![&ae, &bb, c, e]),
        (-6, vec![&ae, &bb, &dd]),
        (-80, vec![&ae, &bd, &cc]),
        (18, vec![a, &bd, c, &dd]),
        (16, vec![&ae, &cc, &cc]),
        (-4, vec![a, &cc, c, &dd]),
        (-27, vec![&bb, &bb, &ee]),
        (18, vec![&bb, &bd, c, e]),
        (-4, vec![&bb, &bd, &dd]),
        (-4, vec![&bb, &cc, c, e]),
        (1, vec![&bb, &cc, &dd]),
    ];
    let mut sum: IPoly = Vec::new();
    for (k, fs) in &terms {
        let p = fs[1..].iter().fold(fs[0].clone(), |acc, f| mul(&acc, f));
        if sum.len() < p.len() {
            sum.resize(p.len(), BigInt::from(0));
        }
        for (s, x) in sum.iter_mut().zip(p) {
            *s += x * *k;
        }
    }
    let d6 = den.pow(6);
    trim(sum.into_iter().map(|x| R::new(x, d6.clone())).collect())
}

/// The same in rationals, product by product (the reference of the test
/// below).
#[cfg(test)]
fn quartic_discriminant_rational(h: &[Poly; 5]) -> Poly {
    let (a, b, c, d, e) = (&h[4], &h[3], &h[2], &h[1], &h[0]);
    let term = |k: i64, fs: &[&Poly]| {
        pscale(
            &fs.iter().fold(vec![int(1)], |acc, p| pmul(&acc, p)),
            &int(k),
        )
    };
    let terms = [
        term(256, &[a, a, a, e, e, e]),
        term(-192, &[a, a, b, d, e, e]),
        term(-128, &[a, a, c, c, e, e]),
        term(144, &[a, a, c, d, d, e]),
        term(-27, &[a, a, d, d, d, d]),
        term(144, &[a, b, b, c, e, e]),
        term(-6, &[a, b, b, d, d, e]),
        term(-80, &[a, b, c, c, d, e]),
        term(18, &[a, b, c, d, d, d]),
        term(16, &[a, c, c, c, c, e]),
        term(-4, &[a, c, c, c, d, d]),
        term(-27, &[b, b, b, b, e, e]),
        term(18, &[b, b, b, c, d, e]),
        term(-4, &[b, b, b, d, d, d]),
        term(-4, &[b, b, c, c, c, e]),
        term(1, &[b, b, c, c, d, d]),
    ];
    trim(terms.iter().fold(Vec::new(), |acc, t| padd(&acc, t)))
}

/// The critical values of `u` (the discriminant's distinct real roots) in a
/// chart whose antipode is none; a torus meridian on the quadric refused.
fn critical(m: &Meeting) -> Result<(Chart, Vec<AlgebraicRoot>)> {
    let cv = Chart {
        c0: int(1),
        s0: zero(),
    };
    for k in 0..16i64 {
        let base = if k == 0 {
            [int(1), zero()]
        } else {
            circle_point(&[zero(), zero()], &int(1), &(int(2 * k + 1) / int(37)))
        };
        let cu = Chart {
            c0: base[0].clone(),
            s0: base[1].clone(),
        };
        let h = h_coefficients(&m.g, &cu, &cv);
        let disc = quartic_discriminant(&h);
        if disc.is_empty() {
            // A double root at every u: tangent along a curve.
            return Err(tangent());
        }
        if disc.len() < 25 {
            continue;
        }
        // A whole tube's circle on the quadric (every coefficient zero at
        // one u): a component no line of u crosses.
        let mut common: Option<IntPolynomial> = None;
        for p in h.iter().filter(|p| !p.is_empty()) {
            let ip = IntPolynomial::from_rationals(p);
            common = Some(match common {
                None => ip,
                Some(c) => c.gcd(&ip),
            });
        }
        if let Some(c) = common.filter(|c| !c.is_constant()) {
            let poly: Poly = trim(c.0.iter().map(|x| R::from_integer(x.clone())).collect());
            if let Ok((_, rs)) = super::turned::roots_repeated(&poly) {
                if !rs.is_empty() {
                    return Err(Error::OutOfDomain(
                        "a torus's tube circle on the other input's surface (S9d.4b.2)",
                    ));
                }
            }
        }
        let rs = if disc.len() == 1 {
            Vec::new()
        } else {
            super::turned::roots_repeated(&disc)?
                .1
                .into_iter()
                .map(|(r, _)| r)
                .collect()
        };
        return Ok((cu, rs));
    }
    Err(limit(
        "a torus meeting's critical values at every chart's antipode",
    ))
}

/// Every critical value's points regular: over the line of `u` at it, each
/// box clear of `G`, of `G_v` or of `G_u` (a turning point, not a singular
/// point of the meeting: the surfaces tangent).
fn check_regular(m: &Meeting, chart: &Chart, crit: &mut [AlgebraicRoot]) -> Result<()> {
    let [g, gu, gv] = &m.cert;
    for r in crit.iter_mut() {
        r.refine_for_signs(100);
        let (lo, hi) = r.isolator();
        let a = angle_r(&chart.at(lo));
        let b = near(angle_r(&chart.at(hi)), a);
        let u = [a.min(b) - 1e-12, a.max(b) + 1e-12];
        let ok = boxes_clear(u, [-PI, PI], |u, v| {
            g.clear(u, v) || gv.clear(u, v) || gu.clear(u, v)
        });
        if !ok {
            return Err(tangent());
        }
    }
    Ok(())
}

// ------------------------------------------------------------ lines

/// A line of `u` at a rational direction between critical values: the
/// meeting's roots on it, exactly and in binary64.
struct Line {
    p: [R; 2],
    pf: f64,
    chart: Chart,
    poly: Poly,
    roots: Vec<AlgebraicRoot>,
    v: Vec<f64>,
}

/// A rational strictly between two roots (ascending, disjoint).
fn rational_between_roots(a: &mut AlgebraicRoot, b: &mut AlgebraicRoot) -> Result<R> {
    for steps in [0, 64, 160] {
        if steps > 0 {
            a.refine_for_signs(steps);
            b.refine_for_signs(steps);
        }
        let m = (a.isolator().1 + b.isolator().0) / int(2);
        if a.compare_rational(&m) == Ordering::Less && b.compare_rational(&m) == Ordering::Greater {
            return Ok(m);
        }
    }
    Err(limit(
        "two critical values of a torus meeting within rounding",
    ))
}

fn lines(m: &Meeting, chart: &Chart, crit: &mut [AlgebraicRoot]) -> Result<Vec<Line>> {
    let mut dirs: Vec<[R; 2]> = Vec::new();
    for i in 0..crit.len().saturating_sub(1) {
        let (x, y) = crit.split_at_mut(i + 1);
        let t = rational_between_roots(&mut x[i], &mut y[0])?;
        dirs.push(chart.at(&t));
    }
    // The wrap gap: the chart's antipode (never critical).
    dirs.push([-chart.c0.clone(), -chart.s0.clone()]);
    lines_at(m, dirs)
}

/// The lines of `u` at rational directions.
fn lines_at(m: &Meeting, dirs: Vec<[R; 2]>) -> Result<Vec<Line>> {
    let mut out = Vec::new();
    for p in dirs {
        let (c, poly, rs) = m.roots_at(false, &p)?;
        let v = rs.iter().map(|r| root_angle(&c, r)).collect();
        out.push(Line {
            pf: angle_r(&p),
            p,
            chart: c,
            poly,
            roots: rs,
            v,
        });
    }
    Ok(out)
}

// ------------------------------------------------------------ folds

/// The binary64 value just above pi.
const PI_HI: f64 = 3.1415926535897936;

/// A box of the angles, `[u] x [v]`.
type Bx = [[f64; 2]; 2];

/// The largest side of a turning point's box.
const FOLD: f64 = 1e-6;

/// Where the meeting may turn in `u` (`G = G_v = 0`), S9d.4b.2b's
/// critical values without a discriminant: a subdivision of `[-pi, pi]^2`
/// into boxes clear of `G` or of `G_v` (none there) and boxes under `FOLD`
/// clear of `G_u` (the meeting regular there: a turning point, not a
/// singular point), merged where they touch (across the period too); a box
/// clear of none at `1e-10` a tangency of the surfaces.
fn folds(m: &Meeting) -> Result<Vec<Bx>> {
    let [g, gu, gv] = &m.cert;
    let mut found: Vec<Bx> = Vec::new();
    let mut stack: Vec<Bx> = vec![[[-PI_HI, PI_HI], [-PI_HI, PI_HI]]];
    let mut budget = 200_000usize;
    while let Some(b) = stack.pop() {
        budget = budget
            .checked_sub(1)
            .ok_or(limit("a torus meeting's turning points"))?;
        let (u, v) = (b[0], b[1]);
        if g.clear(u, v) || gv.clear(u, v) {
            continue;
        }
        let (wu, wv) = (u[1] - u[0], v[1] - v[0]);
        if wu.max(wv) <= FOLD && gu.clear(u, v) {
            found.push(b);
            continue;
        }
        if wu.max(wv) < 1e-10 {
            return Err(tangent());
        }
        let r = usize::from(wv > wu);
        let mid = 0.5 * b[r][0] + 0.5 * b[r][1];
        let (mut lo, mut hi) = (b, b);
        lo[r][1] = mid;
        hi[r][0] = mid;
        stack.push(hi);
        stack.push(lo);
    }
    // Clusters of touching boxes, modulo whole turns.
    let mut clusters: Vec<Bx> = Vec::new();
    for b in found {
        let mut cur = b;
        loop {
            let hit = clusters
                .iter()
                .enumerate()
                .find_map(|(i, c)| touching(c, &cur).map(|shift| (i, shift)));
            let Some((i, [su, sv])) = hit else {
                break;
            };
            let c = clusters.swap_remove(i);
            cur = [
                [c[0][0].min(cur[0][0] + su), c[0][1].max(cur[0][1] + su)],
                [c[1][0].min(cur[1][0] + sv), c[1][1].max(cur[1][1] + sv)],
            ];
        }
        clusters.push(cur);
    }
    Ok(clusters)
}

/// The whole turns `[su, sv]` moving `b` to touch `c`, if any.
fn touching(c: &Bx, b: &Bx) -> Option<[f64; 2]> {
    let meets = |x: [f64; 2], y: [f64; 2], s: f64| x[0] <= y[1] + s && y[0] + s <= x[1];
    for su in [0.0, -TAU, TAU] {
        for sv in [0.0, -TAU, TAU] {
            if meets(c[0], b[0], su) && meets(c[1], b[1], sv) {
                return Some([su, sv]);
            }
        }
    }
    None
}

/// Lines of `u` in the gaps between the turning points' boxes (all of
/// them, at their middles), or at `u = 0` without any.
fn fold_lines(m: &Meeting, clusters: &[Bx]) -> Result<Vec<Line>> {
    let margin = 1e-9;
    // The boxes' ranges of u on the circle, from their starts in [-pi, pi).
    let mut spans: Vec<[f64; 2]> = clusters
        .iter()
        .map(|c| {
            let a = (c[0][0] - margin + PI).rem_euclid(TAU) - PI;
            [a, a + (c[0][1] - c[0][0]) + 2.0 * margin]
        })
        .collect();
    if spans.is_empty() {
        return lines_at(m, vec![[int(1), zero()]]);
    }
    spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
    // Merged runs, the last one carried a turn on when it overlaps the
    // first.
    let mut runs: Vec<[f64; 2]> = Vec::new();
    for s in spans {
        match runs.last_mut() {
            Some(r) if s[0] <= r[1] => r[1] = r[1].max(s[1]),
            _ => runs.push(s),
        }
    }
    while runs.len() > 1 {
        let (first, last) = (runs[0], runs[runs.len() - 1]);
        if last[1] >= first[0] + TAU {
            runs.pop();
            runs[0] = [last[0] - TAU, first[1].max(last[1] - TAU)];
        } else {
            break;
        }
    }
    let n = runs.len();
    let mut dirs = Vec::new();
    for i in 0..n {
        let a = runs[i][1];
        let b = runs[(i + 1) % n][0] + if i + 1 == n { TAU } else { 0.0 };
        if b - a > 1e-7 {
            dirs.push(dir_near(0.5 * a + 0.5 * b));
        }
    }
    if dirs.is_empty() {
        return Err(limit("a torus meeting's turning points all round"));
    }
    lines_at(m, dirs)
}

/// Whether a turning point's box lies on a verified piece over `v`: within
/// its window of `u` and its range of `v` by a margin far above their
/// binary64 views' error (the piece's one root there is then the box's).
fn fold_covered(c: &Bx, pieces: &[ToricCrv]) -> bool {
    let margin = 1e-9;
    let inside = |x: [f64; 2], w: [f64; 2]| {
        let (mid, half) = (0.5 * x[0] + 0.5 * x[1], 0.5 * (x[1] - x[0]));
        let mid = near(mid, 0.5 * w[0] + 0.5 * w[1]);
        w[0] < mid - half - margin && mid + half + margin < w[1]
    };
    pieces
        .iter()
        .any(|p| p.over_v && inside(c[0], p.window_f) && p.range_f.is_none_or(|r| inside(c[1], r)))
}

// ------------------------------------------------------------ tracing

/// A traced component: its points `(u, v)` unwrapped, the last the first
/// moved by its winding's whole turns.
struct Comp {
    pts: Vec<[f64; 2]>,
    winding: [i64; 2],
}

impl Comp {
    fn segments(&self) -> usize {
        self.pts.len() - 1
    }

    /// The point at a position along it (any real: whole laps add the
    /// winding).
    fn at(&self, pos: f64) -> [f64; 2] {
        let n = self.segments() as f64;
        let lap = (pos / n).floor();
        let r = pos - lap * n;
        let i = (r.floor() as usize).min(self.segments() - 1);
        let f = r - i as f64;
        let (a, b) = (self.pts[i], self.pts[i + 1]);
        [0, 1].map(|k| a[k] + f * (b[k] - a[k]) + lap * TAU * self.winding[k] as f64)
    }
}

/// Traces the component through root `ri` of line `li`, marking every line
/// root it passes (`seen`), until it returns there.
fn trace(
    m: &Meeting,
    lines: &[Line],
    seen: &mut [Vec<bool>],
    li: usize,
    ri: usize,
) -> Result<Comp> {
    let start = [lines[li].pf, lines[li].v[ri]];
    let dir_at = |p: [f64; 2]| -> Option<[f64; 2]> {
        let [_, gu, gv] = m.eval(p[0], p[1]);
        let n = gu.hypot(gv);
        (n > 0.0 && n.is_finite()).then(|| [gv / n, -gu / n])
    };
    let mut p = start;
    let mut d = dir_at(p).ok_or(tangent())?;
    let mut h = 1e-2;
    let mut pts = vec![start];
    seen[li][ri] = true;
    for _ in 0..400_000 {
        let pred = [p[0] + h * d[0], p[1] + h * d[1]];
        let mut q = pred;
        let mut converged = false;
        for _ in 0..8 {
            let [g, gu, gv] = m.eval(q[0], q[1]);
            let n2 = gu * gu + gv * gv;
            if n2 <= 0.0 || n2.is_nan() {
                break;
            }
            q = [q[0] - g * gu / n2, q[1] - g * gv / n2];
            if g.abs() / n2.sqrt() < 1e-12 {
                converged = true;
                break;
            }
        }
        let moved = (q[0] - pred[0]).hypot(q[1] - pred[1]);
        let dn = dir_at(q);
        let good =
            converged && moved <= 0.1 * h && dn.is_some_and(|e| e[0] * d[0] + e[1] * d[1] > 0.995);
        if !good {
            h *= 0.5;
            if h < 1e-10 {
                return Err(limit("a torus meeting's trace stalls"));
            }
            continue;
        }
        let dn = dn.expect("a direction");
        // Crossings of the lines of u in (p, q].
        let (u0, u1) = (p[0], q[0]);
        let (lo, hi) = (u0.min(u1), u0.max(u1));
        for (l, line) in lines.iter().enumerate() {
            let kmin = ((lo - line.pf) / TAU).ceil() as i64;
            let kmax = ((hi - line.pf) / TAU).floor() as i64;
            for k in kmin..=kmax {
                let uc = line.pf + TAU * k as f64;
                if uc == u0 || uc < lo || uc > hi {
                    continue;
                }
                let f = (uc - u0) / (u1 - u0);
                let mut vc = p[1] + f * (q[1] - p[1]);
                for _ in 0..4 {
                    let [g, _, gv] = m.eval(uc, vc);
                    if gv != 0.0 {
                        vc -= g / gv;
                    }
                }
                let hit = line
                    .v
                    .iter()
                    .position(|v| (near(*v, vc) - vc).abs() < 1e-7)
                    .ok_or(limit("a torus meeting's trace off its roots"))?;
                if l == li && hit == ri {
                    let end = [uc, near(line.v[ri], vc)];
                    pts.push(end);
                    let winding = [
                        ((end[0] - start[0]) / TAU).round() as i64,
                        ((end[1] - start[1]) / TAU).round() as i64,
                    ];
                    // Exactly the start, moved by whole turns.
                    let n = pts.len() - 1;
                    pts[n] = [
                        start[0] + TAU * winding[0] as f64,
                        start[1] + TAU * winding[1] as f64,
                    ];
                    if n < 3 {
                        return Err(limit("a torus meeting's trace closes at once"));
                    }
                    return Ok(Comp { pts, winding });
                }
                if seen[l][hit] {
                    return Err(limit("a torus meeting's trace crosses another"));
                }
                seen[l][hit] = true;
            }
        }
        pts.push(q);
        p = q;
        d = dn;
        h = (h * 1.5).min(0.05);
    }
    Err(limit("a torus meeting's trace never closes"))
}

// ------------------------------------------------------------ planning

/// A run of a component between positions `a < b`, a graph over `u` or
/// over `v`; a ring closes on itself.
#[derive(Debug, Clone, Copy)]
struct Run {
    a: f64,
    b: f64,
    over_v: bool,
    ring: bool,
}

/// The parameter (`u`, or `v` when `over_v`) and the other angle at a
/// position, the other angle settled on the meeting by Newton's steps.
fn param_at(m: &Meeting, comp: &Comp, over_v: bool, pos: f64) -> (f64, f64) {
    let [u, v] = comp.at(pos);
    let (p, mut s) = if over_v { (v, u) } else { (u, v) };
    for _ in 0..6 {
        let [g, gu, gv] = if over_v { m.eval(s, p) } else { m.eval(p, s) };
        let d = if over_v { gu } else { gv };
        if d == 0.0 {
            break;
        }
        let step = g / d;
        if !step.is_finite() || step.abs() > 0.1 {
            break;
        }
        s -= step;
    }
    (p, s)
}

/// The component's runs: graphs over `u` where its slope in the angles is
/// at most one (`|G_u| <= |G_v|`: its turning points in `v`, `G_u = 0`,
/// among them), over `v` elsewhere (those in `u`), switched where the
/// slope is one, each graph's series as far from its own turning points as
/// the other's; runs over two radians of their parameter split.
fn plan(m: &Meeting, comp: &Comp) -> Vec<Run> {
    let n = comp.segments();
    let slope = |p: [f64; 2]| {
        let [_, gu, gv] = m.eval(p[0], p[1]);
        gv.abs() - gu.abs()
    };
    let d: Vec<f64> = comp.pts.iter().map(|p| slope(*p)).collect();
    let mut switches: Vec<f64> = Vec::new();
    for i in 0..n {
        let (a, b) = (d[i], d[i + 1]);
        if (a >= 0.0) != (b >= 0.0) {
            switches.push(i as f64 + a / (a - b));
        }
    }
    let nf = n as f64;
    let mut runs: Vec<Run> = if switches.is_empty() {
        // One kind all round: a ring over its parameter when it winds once
        // in it and not in the other, else three runs.
        let over_v = d[0] < 0.0;
        let ring = if over_v {
            comp.winding[0] == 0 && comp.winding[1].abs() == 1
        } else {
            comp.winding[1] == 0 && comp.winding[0].abs() == 1
        };
        if ring {
            vec![Run {
                a: 0.0,
                b: nf,
                over_v,
                ring: true,
            }]
        } else {
            (0..3)
                .map(|i| Run {
                    a: nf * i as f64 / 3.0,
                    b: nf * (i + 1) as f64 / 3.0,
                    over_v,
                    ring: false,
                })
                .collect()
        }
    } else {
        let s = switches.len();
        (0..s)
            .map(|j| {
                let a = switches[j];
                let b = if j + 1 == s {
                    switches[0] + nf
                } else {
                    switches[j + 1]
                };
                Run {
                    a,
                    b,
                    over_v: slope(comp.at(0.5 * a + 0.5 * b)) < 0.0,
                    ring: false,
                }
            })
            .collect()
    };
    // Runs over two radians of their parameter split evenly.
    let mut out = Vec::new();
    for r in runs.drain(..) {
        if r.ring {
            out.push(r);
            continue;
        }
        let span = (param_at(m, comp, r.over_v, r.b).0 - param_at(m, comp, r.over_v, r.a).0).abs();
        let parts = (span / 2.0).ceil().max(1.0) as usize;
        for i in 0..parts {
            out.push(Run {
                a: r.a + (r.b - r.a) * i as f64 / parts as f64,
                b: r.a + (r.b - r.a) * (i + 1) as f64 / parts as f64,
                ..r
            });
        }
    }
    out
}

/// Every root's binary64 angle at a rational parameter (repeated roots
/// once).
fn all_roots(m: &Meeting, over_v: bool, p: &[R; 2]) -> Option<Vec<f64>> {
    let form = m.g.at(over_v, p);
    let chart = clear_chart(&form).ok()?;
    let poly = trim(form.poly(&chart));
    if poly.len() < 2 {
        return Some(Vec::new());
    }
    let (_, rs) = super::turned::roots_repeated(&poly).ok()?;
    Some(rs.iter().map(|(r, _)| root_angle(&chart, r)).collect())
}

/// A run's window of the other angle (binary64, ascending, under a turn):
/// between its own values and the other roots' at samples along it.
fn window(m: &Meeting, comp: &Comp, run: &Run) -> Option<[f64; 2]> {
    let k = if run.ring { 48 } else { 16 };
    let (mut prev, mut next) = (f64::NEG_INFINITY, f64::INFINITY);
    let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..=k {
        let pos = run.a + (run.b - run.a) * i as f64 / k as f64;
        let (p, s) = param_at(m, comp, run.over_v, pos);
        let xs = all_roots(m, run.over_v, &dir_near(p))?;
        let own = xs.iter().position(|x| (near(*x, s) - s).abs() < 1e-6)?;
        let (mut pv, mut nx) = (s - TAU, s + TAU);
        for (j, x) in xs.iter().enumerate() {
            if j == own {
                continue;
            }
            let d = (x - s).rem_euclid(TAU);
            nx = nx.min(s + d);
            pv = pv.max(s + d - TAU);
        }
        prev = prev.max(pv);
        next = next.min(nx);
        least = least.min(s);
        most = most.max(s);
    }
    if !(prev < least && most < next) {
        return None;
    }
    let (mut w0, mut w1) = (0.5 * (prev + least), 0.5 * (most + next));
    let cap = TAU - 0.25;
    if w1 - w0 > cap {
        let spare = (cap - (most - least)) / 2.0;
        if spare <= 0.0 {
            return None;
        }
        w0 = w0.max(least - spare);
        w1 = w1.min(most + spare);
    }
    Some([w0, w1])
}

/// A component's pieces, or the run to split.
enum Built {
    Done(Vec<ToricCrv>, Vec<QV>),
    Split(usize),
}

fn build(k: usize, m: &Arc<Meeting>, comp: &Comp, runs: &[Run]) -> Result<Built> {
    let mut windows = Vec::new();
    for (i, r) in runs.iter().enumerate() {
        match window(m, comp, r) {
            Some(w) => windows.push(w),
            None => {
                return Ok(Built::Split(i));
            }
        }
    }
    if runs.len() == 1 && runs[0].ring {
        let piece = ToricCrv::new(k, m, runs[0].over_v, windows[0], None)?;
        return Ok(if piece.verified() {
            Built::Done(vec![piece], Vec::new())
        } else {
            Built::Split(0)
        });
    }
    let r = runs.len();
    // Switches: at each run's start, on the graph over u beside it where
    // there is one.
    let mut switches: Vec<MeetPoint> = Vec::new();
    for j in 0..r {
        let before = (j + r - 1) % r;
        let d = if !runs[before].over_v {
            before
        } else if !runs[j].over_v {
            j
        } else {
            before
        };
        let (p, _) = param_at(m, comp, runs[d].over_v, runs[j].a);
        let probe = ToricCrv::new(k, m, runs[d].over_v, windows[d], None)?;
        match probe.meet_at(&dir_near(p)) {
            Some(x) => switches.push(x),
            None => {
                return Ok(Built::Split(d));
            }
        }
    }
    let mut pieces = Vec::new();
    for j in 0..r {
        let run = &runs[j];
        let (s0, s1) = (&switches[j], &switches[(j + 1) % r]);
        let dir = |x: &MeetPoint| if run.over_v { x.v.clone() } else { x.u.clone() };
        let (a, b) = (dir(s0), dir(s1));
        let rising =
            param_at(m, comp, run.over_v, run.b).0 > param_at(m, comp, run.over_v, run.a).0;
        let range = if rising { [a, b] } else { [b, a] };
        let piece = ToricCrv::new(k, m, run.over_v, windows[j], Some(range))?;
        if !piece.holds_point(s0) || !piece.holds_point(s1) || !piece.verified() {
            return Ok(Built::Split(j));
        }
        pieces.push(piece);
    }
    Ok(Built::Done(
        pieces,
        switches.into_iter().map(|x| x.x).collect(),
    ))
}

/// A component's verified pieces, runs split where one fails.
fn pieces_of(k: usize, m: &Arc<Meeting>, comp: &Comp) -> Result<(Vec<ToricCrv>, Vec<QV>)> {
    let mut runs = plan(m, comp);
    for _ in 0..48 {
        match build(k, m, comp, &runs)? {
            Built::Done(p, s) => return Ok((p, s)),
            Built::Split(i) => {
                let r = runs[i];
                if r.b - r.a < 1e-6 {
                    break;
                }
                let mid = 0.5 * r.a + 0.5 * r.b;
                let halves = [
                    Run {
                        b: mid,
                        ring: false,
                        ..r
                    },
                    Run {
                        a: mid,
                        ring: false,
                        ..r
                    },
                ];
                runs.splice(i..=i, halves);
            }
        }
    }
    Err(limit("a torus meeting's piece left unverified"))
}

// ------------------------------------------------------------ the pair

/// How a whole torus (operand `k`'s model `t`) meets a quadric face of the
/// other input: rings or pieces with their switches, or apart.
pub(super) fn torus_quadric(k: usize, t: &Prism, other: &Other) -> Result<CylPair> {
    torus_far(k, t, &Far::Quadric(Box::new(other.clone())))
}

/// How two whole tori meet (S9d.4b.2b): operand `k`'s model `t` the
/// carrier, `o` the other's.
pub(super) fn torus_torus(k: usize, t: &Prism, o: &Prism) -> Result<CylPair> {
    let ring = o.ring.as_ref().expect("a torus");
    // Equal radii, centres within rounding and axes within rounding of
    // parallel (a frame's normal normalized again, or one normal with the
    // axes turned, rounded): one surface within rounding, whose meeting is
    // no binary64 curve (its projections were left unpinned), as one
    // surface exactly is.
    let own = t.ring.as_ref().expect("a torus");
    let d = sub(&t.f.o, &o.f.o);
    let size = &own.big + &own.small;
    if own.big == ring.big
        && own.small == ring.small
        && dot(&d, &d) * int(10).pow(24) <= &size * &size
        && super::meet::within_rounding_of_parallel(&t.f.n, &o.f.n)
    {
        return Err(Error::Degenerate("two tori within rounding of one surface"));
    }
    let far = Far::Torus {
        f: Box::new(o.f.clone()),
        big: ring.big.clone(),
        small: ring.small.clone(),
    };
    torus_far(k, t, &far)
}

fn torus_far(k: usize, t: &Prism, other: &Far) -> Result<CylPair> {
    let ring = t.ring.as_ref().expect("a torus");
    let m = Arc::new(Meeting::new(&t.f, &ring.big, &ring.small, other));
    if m.g.is_zero() {
        return Err(tangent());
    }
    let mixed = |pieces: Vec<ToricCrv>, switches: Vec<QV>| {
        if pieces.is_empty() {
            return CylPair::Apart;
        }
        CylPair::Mixed(Box::new(super::spheres::Mixed {
            pieces: pieces
                .into_iter()
                .map(|p| Crv::Toric(Box::new(p)))
                .collect(),
            switches,
        }))
    };
    if !m.g.depends_on_u() {
        return Ok(mixed(coaxial(k, &m)?, Vec::new()));
    }
    // Of degree two in each angle, the critical values exactly (the
    // discriminant); two tori's in frames not exactly orthonormal, of
    // degree four, by subdivision (their turning points' boxes).
    let (lines, turns) = if m.g.degree(false) <= 2 && m.g.degree(true) <= 2 {
        let (chart, mut crit) = critical(&m)?;
        check_regular(&m, &chart, &mut crit)?;
        (lines(&m, &chart, &mut crit)?, Vec::new())
    } else {
        let turns = folds(&m)?;
        (fold_lines(&m, &turns)?, turns)
    };
    let mut seen: Vec<Vec<bool>> = lines.iter().map(|l| vec![false; l.v.len()]).collect();
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for li in 0..lines.len() {
        for ri in 0..lines[li].v.len() {
            if seen[li][ri] {
                continue;
            }
            let comp = trace(&m, &lines, &mut seen, li, ri)?;
            let (p, s) = pieces_of(k, &m, &comp)?;
            pieces.extend(p);
            switches.extend(s);
        }
    }
    // Every root on every line lies on a verified piece.
    for line in &lines {
        for root in &line.roots {
            let x = m.exact_point(false, &line.p, &line.chart, &line.poly, root.clone())?;
            if !pieces.iter().any(|p| p.holds_point(&x)) {
                return Err(limit("a torus meeting not covered by its pieces"));
            }
        }
    }
    // Every turning point on a verified piece: a component no line crosses
    // turns in u, so none is missed.
    if !turns.iter().all(|c| fold_covered(c, &pieces)) {
        return Err(limit("a torus meeting's turning point not covered"));
    }
    Ok(mixed(pieces, switches))
}

/// A coaxial quadric's rings: circles at the roots in `v`, rings over `u`.
fn coaxial(k: usize, m: &Arc<Meeting>) -> Result<Vec<ToricCrv>> {
    let (chart, _, rs) = m.roots_at(false, &[int(1), zero()])?;
    let mut vs: Vec<f64> = rs.iter().map(|r| root_angle(&chart, r)).collect();
    vs.sort_by(f64::total_cmp);
    let n = vs.len();
    let mut out = Vec::new();
    for i in 0..n {
        let v = vs[i];
        let window = if n == 1 {
            [v - 3.0, v + 3.0]
        } else {
            let prev = v - (v - vs[(i + n - 1) % n]).rem_euclid(TAU);
            let next = v + (vs[(i + 1) % n] - v).rem_euclid(TAU);
            [0.5 * (prev + v), 0.5 * (v + next)]
        };
        let piece = ToricCrv::new(k, m, false, window, None)?;
        if !piece.verified() {
            return Err(limit("a torus's coaxial circle left unverified"));
        }
        out.push(piece);
    }
    Ok(out)
}

// ------------------------------------------------------------ edges

/// Where the circle `c + a cos + b sin` meets a whole torus (model `t`):
/// the roots of the torus's quartic along it, a polynomial of degree eight
/// in the half-angle tangent (`Q(alpha)`); a tangency refused.
pub(super) fn conic_torus(c: &V, a: &V, b: &V, t: &Prism) -> Result<EdgeMeet> {
    let ring = t.ring.as_ref().expect("a torus");
    let f = &t.f;
    let (l0, la, lb) = (f.local(c), f.local_dir(a), f.local_dir(b));
    let lin: Vec<super::turned::Lin> = (0..3)
        .map(|i| [l0[i].clone(), la[i].clone(), lb[i].clone()])
        .collect();
    let kk = &ring.big * &ring.big - &ring.small * &ring.small;
    let mut s = square_sum(&lin);
    s.add_const(&kk);
    let four = int(4) * &ring.big * &ring.big;
    let form = s.mul(&s).sub(&square_sum(&lin[..2]).scaled(&four));
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let poly = trim(form.poly(&chart));
    if poly.is_empty() && form.value(&[int(-1), zero()]) == zero() {
        return Ok(EdgeMeet::Along);
    }
    let point = |cs: &[Qd; 2]| -> QV {
        [0, 1, 2].map(|j| {
            Qd::rat(c[j].clone())
                .add(&cs[0].scale(&a[j]))
                .add(&cs[1].scale(&b[j]))
        })
    };
    let mut out = Vec::new();
    if form.value(&[int(-1), zero()]) == zero() {
        let cs = [Qd::rat(int(-1)), Qd::rat(zero())];
        out.push((Pos::Ang(cs.clone()), point(&cs)));
    }
    if poly.len() >= 2 {
        for root in roots(&poly)? {
            let g = Arc::new(Gen::new(poly.clone(), root));
            let t = K::generator(&g);
            let den = t.mul(&t).add(&K::Rat(int(1)));
            let inv = den
                .recip()
                .ok_or(limit("a circle's crossing at infinity"))?;
            let cs = [
                Qd::of(K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv)),
                Qd::of(t.scale(&int(2)).mul(&inv)),
            ];
            out.push((Pos::Ang(cs.clone()), point(&cs)));
        }
    }
    Ok(EdgeMeet::Points(out))
}

/// Where a sphere's own circle meets a whole torus: on a basis of equal
/// lengths whose scale to the radius is rational (a whole sphere's great
/// circle, S9d.2) as a conic; any other (a surd radius, unequal axes) by
/// S9d.4c's resultant, a crossing within the resolution `res` of a tangency
/// `Degenerate`.
pub(super) fn circ_torus(circ: &super::sphere::Circ, t: &Prism, res: f64) -> Result<EdgeMeet> {
    let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
    let Some(s) = (xx == yy && dot(&circ.x, &circ.y) == zero())
        .then(|| rational_sqrt(&(&circ.r2 / &xx)))
        .flatten()
    else {
        return super::torus_parts::circ_torus(circ, t, res);
    };
    let (a, b) = (scale(&circ.x, &s), scale(&circ.y, &s));
    match conic_torus(&circ.c, &a, &b, t)? {
        EdgeMeet::Points(p) => Ok(EdgeMeet::Points(
            p.into_iter()
                .map(|(pos, x)| {
                    let Pos::Ang(cs) = pos else {
                        unreachable!("a conic's place")
                    };
                    (Pos::Ang([cs[0].scale(&s), cs[1].scale(&s)]), x)
                })
                .collect(),
        )),
        other => Ok(other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_discriminant_is_the_rational_one() {
        // Coefficient polynomials of mixed degrees and denominators (a
        // small linear congruential stream), zero ones among them.
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) as i64
        };
        for case in 0..40 {
            let h: [Poly; 5] = std::array::from_fn(|j| {
                if (case + j) % 7 == 3 {
                    return Vec::new();
                }
                let len = 1 + (next() % 5) as usize;
                trim(
                    (0..len)
                        .map(|_| int(next() % 2001 - 1000) / int(1 + next() % 97))
                        .collect(),
                )
            });
            assert_eq!(quartic_discriminant(&h), quartic_discriminant_rational(&h));
        }
    }
}
