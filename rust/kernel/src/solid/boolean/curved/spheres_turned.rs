//! S9d.2c: a sphere against a cylinder where a frame is turned
//! (REVIEW_NOTES.md, S9d.2c refined). A cap's or zone's own circles on a
//! turned frame's stored axis lie on bases of unequal lengths: the circle
//! `c + dx x + dy y`, `xx dx^2 + 2 xy dx dy + yy dy^2 = r2`, meets a quadric
//! where, in coordinates `(s, t)` turned from `(dx, dy)` by a rational
//! rotation, the resultant in `t` of the two quadratics vanishes (a quartic
//! in `s`), each meeting's `t` their common root: a point of `Q(alpha)`.
//! A turned cylinder's circle at height `w` is an ellipse in the world, so
//! the sphere's function on it, `F(u, w)`, is of degree two in `(cos u, sin
//! u)`: a loop's piece over the height is the one root of `F(., w)` in its
//! branch's half-turn, verified exactly by the half-turn's boundary quartic
//! `E(w)`, the discriminant in `w` of `F`'s quartic in the half-angle
//! tangent, and a root count.
use super::meet::{tangency, EdgeMeet, Pos};
use super::num::*;
use super::procedural::Other;
use super::sphere::Circ;
use super::spheres::RiseCrv;
use super::turned::{padd, pmul, pscale, roots, roots_repeated, trim, Chart, Form, Poly};
use crate::polynomial::real::{AlgebraicRoot, IntPolynomial};
use crate::solid::split::{q, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::Arc;

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

// ------------------------------------------------------------ circles

/// A quadratic `c0 + c1 X + c2 Y + c3 X^2 + c4 X Y + c5 Y^2`.
type Quad2 = [R; 6];

/// The quadratic in coordinates `(s, t)` turned by `(cr, sr)`: `X = cr s -
/// sr t`, `Y = sr s + cr t`; as `a2 t^2 + a1(s) t + a0(s)`.
fn in_t(c: &Quad2, cr: &R, sr: &R) -> [Poly; 3] {
    let (cc, ss, cs) = (cr * cr, sr * sr, cr * sr);
    let s1 = &c[1] * cr + &c[2] * sr;
    let t1 = &c[2] * cr - &c[1] * sr;
    let s2 = &c[3] * &cc + &c[4] * &cs + &c[5] * &ss;
    let st = int(-2) * &c[3] * &cs + &c[4] * (&cc - &ss) + int(2) * &c[5] * &cs;
    let t2 = &c[3] * &ss - &c[4] * &cs + &c[5] * &cc;
    [
        trim(vec![c[0].clone(), s1, s2]),
        trim(vec![t1, st]),
        trim(vec![t2]),
    ]
}

/// The resultant in `t` of two quadratics `a2 t^2 + a1 t + a0` (`a2` a
/// nonzero constant), a polynomial in `s`.
fn resultant(a: &[Poly; 3], b: &[Poly; 3]) -> Poly {
    let x = padd(&pmul(&a[2], &b[0]), &pscale(&pmul(&a[0], &b[2]), &int(-1)));
    let y = padd(&pmul(&a[2], &b[1]), &pscale(&pmul(&a[1], &b[2]), &int(-1)));
    let z = padd(&pmul(&a[1], &b[0]), &pscale(&pmul(&a[0], &b[1]), &int(-1)));
    trim(padd(&pmul(&x, &x), &pscale(&pmul(&y, &z), &int(-1))))
}

/// A polynomial's value at a field's number.
fn peval_k(p: &Poly, x: &K) -> K {
    p.iter()
        .rev()
        .fold(K::Rat(zero()), |acc, c| acc.mul(x).add(&K::Rat(c.clone())))
}

/// The rotations tried in turn (rational points of the unit circle).
fn rotations() -> [[R; 2]; 4] {
    let r = |a: i64, b: i64, c: i64| [R::new(a.into(), c.into()), R::new(b.into(), c.into())];
    [r(1, 0, 1), r(3, 4, 5), r(5, 12, 13), r(8, 15, 17)]
}

/// The count of distinct real roots, `None` where one is repeated.
fn count(p: &Poly) -> Result<Option<usize>> {
    match roots(p) {
        Ok(rs) => Ok(Some(rs.len())),
        Err(Error::Degenerate(_)) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Where a sphere's circle of unequal axes (a turned cap's rim or split)
/// meets a quadric (S9d.2c): each crossing's `(dx, dy)` place and point in
/// `Q(alpha)`. A tangency, or a crossing within the resolution `res` of
/// one, is `Degenerate`.
pub(super) fn circ_ellipse(circ: &Circ, o: &Other, res: f64) -> Result<EdgeMeet> {
    let (xx, yy, xy) = (
        dot(&circ.x, &circ.x),
        dot(&circ.y, &circ.y),
        dot(&circ.x, &circ.y),
    );
    // The circle: xx X^2 + 2 xy X Y + yy Y^2 - r2.
    let ellipse: Quad2 = [-circ.r2.clone(), zero(), zero(), xx, int(2) * xy, yy];
    // The quadric along the circle's plane: sum (a_i + X lx_i + Y ly_i)^2 -
    // (ra + X rx + Y ry)^2.
    let mut quad: Quad2 = std::array::from_fn(|_| zero());
    let mut add_sq = |l: [R; 3], sign: &R| {
        let [a, x, y] = l;
        quad[0] += sign * &a * &a;
        quad[1] += sign * int(2) * &a * &x;
        quad[2] += sign * int(2) * &a * &y;
        quad[3] += sign * &x * &x;
        quad[4] += sign * int(2) * &x * &y;
        quad[5] += sign * &y * &y;
    };
    for (g, e) in o.g.iter().zip(&o.e) {
        add_sq(
            [dot(g, &circ.c) - e, dot(g, &circ.x), dot(g, &circ.y)],
            &int(1),
        );
    }
    let rad = o.radius_lin(&circ.c, &circ.x, &circ.y);
    let rho = if rad[0] < zero() {
        -rad[0].clone()
    } else {
        rad[0].clone()
    };
    add_sq(rad, &int(-1));
    // The band of the quadric's function about zero a displacement of the
    // resolution spans: 2 rho res (its rows about unit).
    let delta = int(2) * (rho + q(res)) * q(res);
    let offset = |k: &R| {
        let mut c = quad.clone();
        c[0] += k;
        c
    };
    let mut tangent = false;
    for [cr, sr] in rotations() {
        let a = in_t(&ellipse, &cr, &sr);
        let b = in_t(&quad, &cr, &sr);
        let res0 = resultant(&a, &b);
        if res0.is_empty() {
            // The circle on the surface.
            return Err(tangency());
        }
        let Some(n) = count(&res0)? else {
            tangent = true;
            continue;
        };
        // An extremum of the quadric's function along the circle within
        // the band: a crossing within the resolution of a tangency.
        for k in [delta.clone(), -delta.clone()] {
            let rk = resultant(&a, &in_t(&offset(&k), &cr, &sr));
            if count(&rk)? != Some(n) {
                return Err(Error::Degenerate(
                    "a sphere's circle within the resolution of tangency to a surface (S9d.2c)",
                ));
            }
        }
        let mut out = Vec::new();
        let mut ok = true;
        for root in roots(&res0)? {
            let g = Arc::new(Gen::new(res0.clone(), root));
            let s = K::generator(&g);
            let ev = |p: &Poly| peval_k(p, &s);
            // b2 A - a2 B = (b2 a1 - a2 b1) t + (b2 a0 - a2 b0).
            let num = ev(&a[2]).mul(&ev(&b[0])).sub(&ev(&a[0]).mul(&ev(&b[2])));
            let den = ev(&a[1]).mul(&ev(&b[2])).sub(&ev(&a[2]).mul(&ev(&b[1])));
            let Some(inv) = den.recip() else {
                ok = false;
                break;
            };
            let t = num.mul(&inv);
            let dx = s.scale(&cr).sub(&t.scale(&sr));
            let dy = s.scale(&sr).add(&t.scale(&cr));
            let p: QV = [0, 1, 2].map(|j| {
                Qd::of(
                    K::Rat(circ.c[j].clone())
                        .add(&dx.scale(&circ.x[j]))
                        .add(&dy.scale(&circ.y[j])),
                )
            });
            out.push((Pos::Ang([Qd::of(dx), Qd::of(dy)]), p));
        }
        if ok {
            return Ok(EdgeMeet::Points(out));
        }
        tangent = true;
    }
    debug_assert!(tangent);
    Err(tangency())
}

// ------------------------------------------------------------ loops

/// `F(cos, sin; w)` of a sphere (`c`, `rr`) on a cylinder's circle `o + r
/// (cos x + sin y) + w n` at height `w`: its terms by exponents of `(cos,
/// sin)`, coefficients polynomials in `w`.
fn f_terms(o: &V, x: &V, y: &V, n: &V, r: &R, c: &V, rr: &R) -> Vec<((u32, u32), Poly)> {
    let d = sub(o, c);
    vec![
        (
            (0, 0),
            trim(vec![dot(&d, &d) - rr * rr, int(2) * dot(&d, n), dot(n, n)]),
        ),
        ((2, 0), trim(vec![r * r * dot(x, x)])),
        ((1, 1), trim(vec![int(2) * r * r * dot(x, y)])),
        ((0, 2), trim(vec![r * r * dot(y, y)])),
        (
            (1, 0),
            trim(vec![int(2) * r * dot(x, &d), int(2) * r * dot(x, n)]),
        ),
        (
            (0, 1),
            trim(vec![int(2) * r * dot(y, &d), int(2) * r * dot(y, n)]),
        ),
    ]
}

fn pvalue(p: &Poly, w: &R) -> R {
    p.iter().rev().fold(zero(), |acc, c| acc * w + c)
}

/// `F(., w)` at a rational height, a quadratic form in `(cos, sin)`.
fn form_at(terms: &[((u32, u32), Poly)], w: &R) -> Form {
    let map: BTreeMap<(u32, u32), R> = terms.iter().map(|(e, p)| (*e, pvalue(p, w))).collect();
    Form::from_terms(map, 2)
}

/// The roots of `F(., w)` on the side `plus` of `(alpha, beta)` at a
/// rational height, as `(cos, sin)` in `Q(alpha)` (a repeated real root, a
/// tangency, refused).
fn half_turn(terms: &[((u32, u32), Poly)], w: &R, ab: &[R; 2], plus: bool) -> Result<Vec<[K; 2]>> {
    let form = form_at(terms, w);
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let p = trim(form.poly(&chart));
    if p.is_empty() {
        return Err(tangency());
    }
    let want = if plus {
        Ordering::Greater
    } else {
        Ordering::Less
    };
    let mut out = Vec::new();
    // The chart's antipode (-1, 0), where the side is beta's sign.
    if form.value(&[int(-1), zero()]) == zero() && sign(&ab[1]) == want {
        out.push([K::Rat(int(-1)), K::Rat(zero())]);
    }
    // alpha sin - beta cos, times 1 + t^2.
    let side = IntPolynomial::from_rationals(&[-ab[1].clone(), int(2) * &ab[0], ab[1].clone()]);
    for root in roots(&p)? {
        if root.sign_polynomial(&side) != want {
            continue;
        }
        let g = Arc::new(Gen::new(p.clone(), root));
        let t = K::generator(&g);
        let inv = t
            .mul(&t)
            .add(&K::Rat(int(1)))
            .recip()
            .ok_or(limit("a height graph's point at infinity"))?;
        out.push([
            K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv),
            t.scale(&int(2)).mul(&inv),
        ]);
    }
    Ok(out)
}

/// A turned piece's point at a rational height: the one root of `F(., w)`
/// on its side (`None` where there is none or more than one).
pub(super) fn rise_at(c: &RiseCrv, w: &R) -> Result<Option<QV>> {
    let terms = f_terms(&c.o, &c.x, &c.y, &c.n, &c.r, &c.c, &c.rr);
    let (alpha, beta, _) = c.coefficients();
    let mut found = half_turn(&terms, w, &[alpha, beta], c.plus)?;
    if found.len() != 1 {
        return Ok(None);
    }
    let [co, si] = found.pop().expect("one root");
    Ok(Some([0, 1, 2].map(|j| {
        Qd::of(
            K::Rat(&c.o[j] + w * &c.n[j])
                .add(&co.scale(&(&c.r * &c.x[j])))
                .add(&si.scale(&(&c.r * &c.y[j]))),
        )
    })))
}

/// A sphere's meeting with a turned cylinder over its height (S9d.2c): `F`'s
/// terms, the half-turns' axis `(alpha, beta)`, and the heights where a
/// root reaches a half-turn's boundary (`E`'s real roots) or two roots meet
/// (the discriminant's: the height graph's turning points).
#[derive(Debug)]
pub(super) struct Height {
    terms: Vec<((u32, u32), Poly)>,
    ab: [R; 2],
    edge: Vec<AlgebraicRoot>,
    turns: Vec<AlgebraicRoot>,
}

impl Height {
    /// The meeting of the sphere (`c`, `rr`) with the cylinder `o + r (cos x
    /// + sin y) + w n`, `(alpha, beta)` the half-turns' axis.
    pub(super) fn new(
        (o, x, y, n): (&V, &V, &V, &V),
        r: &R,
        c: &V,
        rr: &R,
        ab: [R; 2],
    ) -> Result<Self> {
        let terms = f_terms(o, x, y, n, r, c, rr);
        // F (1 + t^2)^2 = sum_j h_j(w) t^j in the chart based at (1, 0).
        let chart = Chart {
            c0: int(1),
            s0: zero(),
        };
        let [cn, sn] = chart.numerators();
        let w2 = vec![int(1), zero(), int(1)];
        let pw = |p: &Poly, k: u32| (0..k).fold(vec![int(1)], |acc, _| pmul(&acc, p));
        let mut h: [Poly; 5] = Default::default();
        for ((i, j), coef) in &terms {
            let tp = pmul(&pmul(&pw(&cn, *i), &pw(&sn, *j)), &pw(&w2, 2 - i - j));
            for (k, tc) in tp.iter().enumerate() {
                h[k] = padd(&h[k], &pscale(coef, tc));
            }
        }
        let disc = super::torus_curved::quartic_discriminant(&h);
        if disc.is_empty() {
            return Err(tangency());
        }
        let turns = if disc.len() == 1 {
            Vec::new()
        } else {
            roots_repeated(&disc)?
                .1
                .into_iter()
                .map(|(r, _)| r)
                .collect()
        };
        // On the boundary directions +-(alpha, beta) / sqrt(rho2): F = P(w) +-
        // l(w) / sqrt(rho2).
        let [alpha, beta] = &ab;
        let rho2 = alpha * alpha + beta * beta;
        let quad = r
            * r
            * (dot(x, x) * alpha * alpha
                + int(2) * dot(x, y) * alpha * beta
                + dot(y, y) * beta * beta)
            / &rho2;
        let d = sub(o, c);
        let p = padd(
            &trim(vec![quad]),
            &trim(vec![dot(&d, &d) - rr * rr, int(2) * dot(&d, n), dot(n, n)]),
        );
        let l = trim(vec![
            int(2) * r * (alpha * dot(x, &d) + beta * dot(y, &d)),
            int(2) * r * (alpha * dot(x, n) + beta * dot(y, n)),
        ]);
        let e = trim(padd(
            &pscale(&pmul(&p, &p), &rho2),
            &pscale(&pmul(&l, &l), &int(-1)),
        ));
        if e.is_empty() {
            return Err(limit("a sphere along a turned cylinder's branch line"));
        }
        let edge = if e.len() == 1 {
            Vec::new()
        } else {
            roots_repeated(&e)?.1.into_iter().map(|(r, _)| r).collect()
        };
        Ok(Self {
            terms,
            ab,
            edge,
            turns,
        })
    }

    /// The height graph's turning points (the ellipse tangent to the
    /// sphere), ascending.
    pub(super) fn turns(&self) -> Vec<AlgebraicRoot> {
        self.turns.clone()
    }

    /// Whether the piece over `lo..=hi` on the side `plus` is verified: no
    /// root on the half-turn's boundary and no double root over its heights,
    /// one root on its side at a rational height inside.
    pub(super) fn verified(&self, lo: &Qd, hi: &Qd, plus: bool) -> Result<bool> {
        let (il, ih) = (lo.interval(), hi.interval());
        let (l, h) = (il.lo().clone(), ih.hi().clone());
        let clear = |rs: &[AlgebraicRoot]| {
            rs.iter().all(|r| {
                r.compare_rational(&l) == Ordering::Less
                    || r.compare_rational(&h) == Ordering::Greater
            })
        };
        if !clear(&self.edge) || !clear(&self.turns) {
            return Ok(false);
        }
        let m = (il.hi() + ih.lo()) / int(2);
        let mq = Qd::rat(m.clone());
        if lo.cmp(&mq) != Ordering::Less || mq.cmp(hi) != Ordering::Less {
            return Ok(false);
        }
        Ok(half_turn(&self.terms, &m, &self.ab, plus)?.len() == 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solid::boolean::curved::model::Affine;
    use crate::solid::boolean::curved::procedural::other_of;
    use crate::{Frame3, Point3, Tolerance, Vec3};

    fn frame(n: [f64; 3], x: [f64; 3]) -> Affine {
        let f = Frame3::new(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            Tolerance::default(),
        )
        .unwrap();
        Affine::new(&f).unwrap()
    }

    fn points(m: &EdgeMeet) -> Vec<[f64; 3]> {
        let EdgeMeet::Points(ps) = m else {
            panic!("points")
        };
        let mut out: Vec<[f64; 3]> = ps.iter().map(|(_, p)| qv_f64(p)).collect();
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out
    }

    /// On a circle of equal axes the resultant finds S9d.2a's surd points.
    #[test]
    fn the_resultant_meets_a_round_circle_as_the_surds_do() {
        // The unit sphere's equator scaled by 2 against the cylinder of
        // radius 1 about x = 1.5 along z: two crossings.
        let circ = Circ {
            c: [zero(), zero(), zero()],
            x: [int(1), zero(), zero()],
            y: [zero(), int(1), zero()],
            r2: int(4),
        };
        let f = frame([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let other = other_of(&f, &[q(1.5), zero()], &int(1));
        let res = 1e-7;
        let a = points(&super::super::spheres::circ_quadric(&circ, &other, res).unwrap());
        let b = points(&circ_ellipse(&circ, &other, res).unwrap());
        assert_eq!(a.len(), 2);
        assert_eq!(a.len(), b.len());
        for (p, r) in a.iter().zip(&b) {
            for k in 0..3 {
                assert!((p[k] - r[k]).abs() < 1e-15, "{p:?} {r:?}");
            }
        }
    }

    /// A circle of unequal axes: its meetings lie on it and on the cylinder
    /// exactly, and a crossing within the resolution of tangency is
    /// refused.
    #[test]
    fn a_circle_of_unequal_axes_meets_a_cylinder_exactly() {
        let circ = Circ {
            c: [zero(), zero(), zero()],
            x: [int(3), zero(), zero()],
            y: [zero(), int(4), zero()],
            r2: int(144),
        };
        // X^2 9 + Y^2 16 = 144: the circle of radius 12 about the origin.
        let f = frame([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let other = other_of(&f, &[int(10), zero()], &int(3));
        let m = circ_ellipse(&circ, &other, 1e-7).unwrap();
        let EdgeMeet::Points(ps) = &m else {
            panic!("points")
        };
        assert_eq!(ps.len(), 2);
        for (_, p) in ps {
            assert!(circ.on(p));
            assert_eq!(other.value(p).sign(), Ordering::Equal);
        }
        // Tangent from outside, 1e-9 short: within the resolution.
        let near = other_of(&f, &[q(15.0 + 1e-9), zero()], &int(3));
        assert!(matches!(
            circ_ellipse(&circ, &near, 1e-7),
            Err(Error::Degenerate(_))
        ));
    }

    /// In an exact frame the discriminant's real roots are S9d.2b's `dw`'s.
    #[test]
    fn an_exact_frames_turning_points_are_the_quartics() {
        let (o, x, y, n) = (
            [zero(), zero(), zero()],
            [int(1), zero(), zero()],
            [zero(), int(1), zero()],
            [zero(), zero(), int(1)],
        );
        let (r, c, rr) = (int(1), [int(-2), q(-0.5), zero()], int(2));
        let d = sub(&o, &c);
        let ab = [int(2) * dot(&x, &d), int(2) * dot(&y, &d)];
        let h = Height::new((&o, &x, &y, &n), &r, &c, &rr, ab.clone()).unwrap();
        // dw = rho2 r^2 - g(w)^2, g = rr^2 - |d|^2 - r^2 - w^2 (n . d = 0).
        let rho2 = &ab[0] * &ab[0] + &ab[1] * &ab[1];
        let g = vec![&rr * &rr - dot(&d, &d) - &r * &r, zero(), int(-1)];
        let dw = padd(&vec![&rho2 * &r * &r], &pscale(&pmul(&g, &g), &int(-1)));
        let want: Vec<f64> = roots(&dw)
            .unwrap()
            .iter()
            .map(|x| crate::solid::split::rational_f64(&super::super::turned::middle(x)))
            .collect();
        let got: Vec<f64> = h
            .turns()
            .iter()
            .map(|x| crate::solid::split::rational_f64(&super::super::turned::middle(x)))
            .collect();
        assert_eq!(want.len(), 2);
        assert_eq!(got.len(), want.len());
        for (a, b) in got.iter().zip(&want) {
            assert!((a - b).abs() < 1e-6, "{got:?} {want:?}");
        }
    }
}
