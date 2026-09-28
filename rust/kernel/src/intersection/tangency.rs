//! Exact tangencies of a torus and a cylinder (S7b.3b of `REVIEW_NOTES.md`).
//!
//! A torus is the pipe of radius `r` about its spine circle and a cylinder the
//! pipe of radius `r_c` about its axis, so they touch exactly where the
//! distance between the spine and the axis has a critical point of value
//! `r + r_c` or `|r - r_c|`: there the segment between the spine point `s` and
//! its foot `l` on the axis is normal to both, and the contact is on it, `r`
//! from `s`. In a rational basis `(U, V)` of the spine's plane,
//! `s = o + u U + v V` satisfies three conics: the circle
//! `|u U + v V|^2 = R^2` (`E0`), the critical-point condition
//! `M (s - o2) . (a x (s - o)) = 0` (`E1`, `M` the projector normal to the
//! axis) and the distance `|M (s - o2)|^2 = k` (`E2`). The critical points
//! are the real roots `u` of the resultant of `E0` and `E1` in `v`, with `v`
//! their common root, rational in `u`; a critical point is a tangency when
//! `E2`'s numerator vanishes there, decided exactly on the algebraic root. The
//! basis is sheared until no critical point has a vanishing denominator.
use super::procedural::{cross, dot, sub, zero, X};
use crate::certified::{Interval as I, Real};
use crate::polynomial::real::{isolate, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::{Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;

/// Univariate rational polynomials, ascending powers.
type P = Vec<R>;

fn rz() -> R {
    R::from_integer(0.into())
}
fn padd(a: &P, b: &P) -> P {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| a.get(i).cloned().unwrap_or_else(rz) + b.get(i).cloned().unwrap_or_else(rz))
        .collect()
}
fn psub(a: &P, b: &P) -> P {
    padd(a, &pscale(b, &R::from_integer((-1).into())))
}
fn pscale(a: &P, k: &R) -> P {
    a.iter().map(|x| x * k).collect()
}
fn pmul(a: &P, b: &P) -> P {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![rz(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}
/// `c0 + c1 u`.
fn lin(c0: R, c1: R) -> P {
    vec![c0, c1]
}
fn enclose(p: &P, u: &I) -> I {
    p.iter().rev().fold(I::exact_f64(0.0), |acc, c| {
        acc.mul(u).add(&I::exact(c.clone()))
    })
}

/// A tangency: the spine point `s`, the contact point, both enclosed.
pub(super) struct Contact {
    pub(super) spine: [I; 3],
    pub(super) point: [I; 3],
}

/// The exact tangencies of the torus (origin `o`, axis `a`, radii `big`,
/// `small`) and the cylinder (origin `o2`, axis `a2`, radius `rc`) whose axes
/// are not coaxial. `ComputationLimit` for an axis meeting the spine with
/// equal radii, or when no shear separates the critical points.
pub(super) fn torus_cylinder(
    o: &X,
    a: &X,
    big: &R,
    small: &R,
    o2: &X,
    a2: &X,
    rc: &R,
) -> Result<Vec<Contact>> {
    let limit = || Error::ComputationLimit("a torus's tangency with a cylinder");
    let u0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        .iter()
        .map(|e| cross(a, &e.map(|k| R::from_integer(k.into()))))
        .find(|u| !u.iter().all(zero))
        .expect("a nonzero axis");
    let v0 = cross(a, &u0);
    let aa2 = dot(a2, a2);
    // M x = x - (x . a2) a2 / |a2|^2, and M x . M y = x . M y.
    let m = |x: &X| -> X {
        let k = dot(x, a2) / &aa2;
        sub(x, &a2.clone().map(|c| c * &k))
    };
    let d = sub(o, o2);
    // The critical values of a contact: (r + r_c)^2, then (r - r_c)^2.
    let ks = {
        let (p, q) = (small + rc, small - rc);
        let mut ks = vec![&p * &p];
        if &q * &q != &p * &p {
            ks.push((small - rc) * (small - rc));
        }
        ks
    };
    for shear in [0i64, 1, -2, 3, -5, 7, 11, -13] {
        // V = v0 + (shear / 7) U: points with one u lie along V, so the
        // shear changes which critical points share a u.
        let sh = R::new(shear.into(), 7.into());
        let u = u0.clone();
        let v = [
            &v0[0] + &u0[0] * &sh,
            &v0[1] + &u0[1] * &sh,
            &v0[2] + &u0[2] * &sh,
        ];
        let (au, av) = (cross(a, &u), cross(a, &v));
        let (mu, mv, md) = (m(&u), m(&v), m(&d));
        // E0: a2 v^2 + a1 v + a0.
        let e0 = [
            lin(-(big * big), rz())
                .into_iter()
                .chain([dot(&u, &u)])
                .collect::<P>(),
            lin(rz(), R::from_integer(2.into()) * dot(&u, &v)),
            vec![dot(&v, &v)],
        ];
        // E1: M(d + u U + v V) . (u A x U + v A x V).
        let e1 = [
            lin(rz(), dot(&md, &au))
                .into_iter()
                .chain([dot(&mu, &au)])
                .collect::<P>(),
            lin(dot(&md, &av), dot(&mu, &av) + dot(&mv, &au)),
            vec![dot(&mv, &av)],
        ];
        // E2 without k: |M (d + u U + v V)|^2.
        let two = R::from_integer(2.into());
        let e2 = [
            vec![dot(&md, &md), &two * dot(&md, &mu), dot(&mu, &mu)],
            lin(&two * dot(&md, &mv), &two * dot(&mu, &mv)),
            vec![dot(&mv, &mv)],
        ];
        let [a0, a1, a2c] = &e0;
        let [b0, b1, b2] = &e1;
        let res = psub(
            &pmul(
                &psub(&pmul(a2c, b0), &pmul(a0, b2)),
                &psub(&pmul(a2c, b0), &pmul(a0, b2)),
            ),
            &pmul(
                &psub(&pmul(a2c, b1), &pmul(a1, b2)),
                &psub(&pmul(a1, b0), &pmul(a0, b1)),
            ),
        );
        let resultant = IntPolynomial::from_rationals(&res);
        if resultant.is_zero() {
            return Err(limit());
        }
        // v = Nv / Dn at a common root.
        let nv = psub(&pmul(a2c, b0), &pmul(a0, b2));
        let dn = psub(&pmul(a1, b2), &pmul(a2c, b1));
        let dn_int = IntPolynomial::from_rationals(&dn);
        let bound = cauchy(&res);
        let mut budget = Budget::new(RootIsolationOptions::default());
        let roots = isolate(&resultant, -bound.clone(), bound, &mut budget)?;
        if roots
            .iter()
            .any(|x| dn_int.is_zero() || x.vanishes_polynomial(&dn_int))
        {
            continue;
        }
        let mut out = Vec::new();
        for mut root in roots {
            let [c0, c1, c2] = &e2;
            for k in &ks {
                // Numerator of E2 - k at v = Nv / Dn.
                let h = padd(
                    &padd(&pmul(c2, &pmul(&nv, &nv)), &pmul(c1, &pmul(&nv, &dn))),
                    &pmul(&psub(c0, &vec![k.clone()]), &pmul(&dn, &dn)),
                );
                let hi = IntPolynomial::from_rationals(&h);
                if !(hi.is_zero() || root.vanishes_polynomial(&hi)) {
                    continue;
                }
                root.refine_for_signs(400);
                let (lo, up) = root.isolator();
                let uu = I::new(lo.clone(), up.clone());
                let vv = enclose(&nv, &uu)
                    .div(&enclose(&dn, &uu))
                    .ok_or_else(limit)?;
                let x = |c: &R| I::exact(c.clone());
                let s: [I; 3] = std::array::from_fn(|i| {
                    x(&o[i]).add(&uu.mul(&x(&u[i]))).add(&vv.mul(&x(&v[i])))
                });
                // The foot l on the axis, and the contact on the segment.
                let w: [I; 3] = std::array::from_fn(|i| s[i].sub(&x(&o2[i])));
                let h = w
                    .iter()
                    .zip(a2)
                    .fold(I::exact_f64(0.0), |acc, (wi, ai)| acc.add(&wi.mul(&x(ai))))
                    .div(&x(&aa2))
                    .ok_or_else(limit)?;
                let gap: [I; 3] = std::array::from_fn(|i| h.mul(&x(&a2[i])).sub(&w[i]));
                let norm = gap
                    .iter()
                    .fold(I::exact_f64(0.0), |acc, g| acc.add(&g.square()))
                    .sqrt();
                // The contact is r from s towards l, except inside a thicker
                // cylinder (k = (r - r_c)^2 with r_c > r), where it is away.
                let away = k != &ks[0] && rc > small;
                let r = if away { x(small).neg() } else { x(small) };
                if norm.sign() != Some(std::cmp::Ordering::Greater) {
                    return Err(limit());
                }
                let point: [I; 3] = std::array::from_fn(|i| {
                    s[i].add(&r.mul(&gap[i]).div(&norm).expect("a positive distance"))
                });
                out.push(Contact { spine: s, point });
            }
        }
        return Ok(out);
    }
    Err(limit())
}

/// A bound on the absolute values of a polynomial's real roots.
fn cauchy(p: &P) -> R {
    let lead = p
        .iter()
        .rev()
        .find(|c| **c != rz())
        .cloned()
        .unwrap_or_else(|| R::from_integer(1.into()));
    let m = p
        .iter()
        .map(|c| {
            let x = c / &lead;
            if x < rz() {
                -x
            } else {
                x
            }
        })
        .fold(rz(), |a, b| if b > a { b } else { a });
    m + R::from_integer(BigInt::from(1))
}

// ------------------------------------------------------------------ tori

/// Polynomials in `(u, v)`: the coefficients of `v^0, v^1, ...`, each a
/// polynomial in `u`.
type B = Vec<P>;

fn bl(c: R, cu: R, cv: R) -> B {
    vec![vec![c, cu], vec![cv]]
}
fn badd(a: &B, b: &B) -> B {
    let n = a.len().max(b.len());
    (0..n)
        .map(|k| {
            padd(
                a.get(k).map_or(&Vec::new(), |x| x),
                b.get(k).map_or(&Vec::new(), |x| x),
            )
        })
        .collect()
}
fn bscale(a: &B, k: &R) -> B {
    a.iter().map(|x| pscale(x, k)).collect()
}
fn bmul(a: &B, b: &B) -> B {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = padd(&out[i + j], &pmul(x, y));
        }
    }
    out
}
fn bdot(a: &[B; 3], b: &[B; 3]) -> B {
    badd(
        &badd(&bmul(&a[0], &b[0]), &bmul(&a[1], &b[1])),
        &bmul(&a[2], &b[2]),
    )
}
/// `b` modulo `a2 v^2 + a1 v + a0` (`a2` a nonzero constant): `c0 + c1 v`.
fn reduce(b: &B, a0: &P, a1: &P, a2: &R) -> (P, P) {
    let mut b = b.clone();
    for n in (2..b.len()).rev() {
        let top = std::mem::take(&mut b[n]);
        let k = pscale(&top, &(R::from_integer((-1).into()) / a2));
        b[n - 1] = padd(&b[n - 1], &pmul(&k, a1));
        b[n - 2] = padd(&b[n - 2], &pmul(&k, a0));
    }
    let get = |k: usize| b.get(k).cloned().unwrap_or_default();
    (get(0), get(1))
}

/// The exact tangencies of two tori (origins `o1`, `o2`, axes `a1`, `a2`,
/// radii `(R1, r1)`, `(R2, r2)`) off a common axis: pairs of spine points at
/// a critical distance `r1 + r2` or `|r1 - r2|`, for the first spine's point
/// `s1 = o1 + u U + v V` on its circle `E0`. With `w = s1 - o2`, `M` the
/// projector normal to `a2` and `Q = k - |w|^2 - R2^2`, the distance to the
/// second spine is `k` and critical along the first at `s1` exactly when
/// `E1 = (w . T1) Q + 2 R2^2 (M w . T1) = 0` and
/// `E2 = Q^2 - 4 R2^2 |M w|^2 = 0` (`T1 = a1 x (s1 - o1)`, and `w . T1` is
/// linear in `u`, `v`). Both are reduced modulo `E0` to `c0 + c1 v` and
/// `d0 + d1 v`; the critical points are the real roots of the resultant
/// `a2 c0^2 - a1 c0 c1 + a0 c1^2`, `v = -c0 / c1`, and a tangency where
/// `d0 c1 - d1 c0` vanishes on the algebraic root. `ComputationLimit` for a
/// critical point on the second torus's axis (`Q = 0`), spines meeting with
/// equal minor radii, tori touching along a curve, or when no shear
/// separates the critical points.
#[allow(clippy::too_many_arguments)]
pub(super) fn torus_torus(
    o1: &X,
    a1: &X,
    big1: &R,
    small1: &R,
    o2: &X,
    a2: &X,
    big2: &R,
    small2: &R,
) -> Result<Vec<Contact>> {
    let limit = || Error::ComputationLimit("the tangency of two tori");
    let u0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        .iter()
        .map(|e| cross(a1, &e.map(|k| R::from_integer(k.into()))))
        .find(|u| !u.iter().all(zero))
        .expect("a nonzero axis");
    let v0 = cross(a1, &u0);
    let aa2 = dot(a2, a2);
    let d = sub(o1, o2);
    let ks = {
        let (p, q) = (small1 + small2, small1 - small2);
        let mut ks = vec![&p * &p];
        if &q * &q != &p * &p {
            ks.push(&q * &q);
        }
        ks
    };
    let two = R::from_integer(2.into());
    let four = R::from_integer(4.into());
    for shear in [0i64, 1, -2, 3, -5, 7, 11, -13] {
        let sh = R::new(shear.into(), 7.into());
        let u = u0.clone();
        let v: X = std::array::from_fn(|i| &v0[i] + &u0[i] * &sh);
        let (au, av) = (cross(a1, &u), cross(a1, &v));
        // E0 = a2 v^2 + a1 v + a0.
        let e0a0: P = vec![-(big1 * big1), rz(), dot(&u, &u)];
        let e0a1: P = vec![rz(), &two * dot(&u, &v)];
        let e0a2 = dot(&v, &v);
        let w: [B; 3] = std::array::from_fn(|i| bl(d[i].clone(), u[i].clone(), v[i].clone()));
        let t1: [B; 3] = std::array::from_fn(|i| bl(rz(), au[i].clone(), av[i].clone()));
        let wa = bdot(&w, &std::array::from_fn(|i| vec![vec![a2[i].clone()]]));
        let mw: [B; 3] = std::array::from_fn(|i| badd(&w[i], &bscale(&wa, &(-(&a2[i] / &aa2)))));
        let ww = bdot(&w, &w);
        let wt = bdot(&w, &t1);
        let mt = bdot(&mw, &t1);
        let mm = bdot(&mw, &mw);
        let r22 = big2 * big2;
        let mut found = Vec::new();
        let mut separated = true;
        for k in &ks {
            let q = badd(
                &vec![vec![k - &r22]],
                &bscale(&ww, &R::from_integer((-1).into())),
            );
            let e1 = badd(&bmul(&wt, &q), &bscale(&mt, &(&two * &r22)));
            let e2 = badd(&bmul(&q, &q), &bscale(&mm, &(-(&four * &r22))));
            let (c0, c1) = reduce(&e1, &e0a0, &e0a1, &e0a2);
            let (d0, d1) = reduce(&e2, &e0a0, &e0a1, &e0a2);
            let (q0, q1) = reduce(&q, &e0a0, &e0a1, &e0a2);
            let res = padd(
                &padd(
                    &pscale(&pmul(&c0, &c0), &e0a2),
                    &pscale(&pmul(&pmul(&e0a1, &c0), &c1), &R::from_integer((-1).into())),
                ),
                &pmul(&e0a0, &pmul(&c1, &c1)),
            );
            let resultant = IntPolynomial::from_rationals(&res);
            if resultant.is_zero() {
                return Err(limit());
            }
            let c1i = IntPolynomial::from_rationals(&c1);
            let mut budget = Budget::new(RootIsolationOptions::default());
            let bound = cauchy(&res);
            let roots = isolate(&resultant, -bound.clone(), bound, &mut budget)?;
            if roots
                .iter()
                .any(|x| c1i.is_zero() || x.vanishes_polynomial(&c1i))
            {
                separated = false;
                break;
            }
            let tangent = IntPolynomial::from_rationals(&psub(&pmul(&d0, &c1), &pmul(&d1, &c0)));
            let on_axis = IntPolynomial::from_rationals(&psub(&pmul(&q0, &c1), &pmul(&q1, &c0)));
            for mut root in roots {
                if !(tangent.is_zero() || root.vanishes_polynomial(&tangent)) {
                    continue;
                }
                if on_axis.is_zero() || root.vanishes_polynomial(&on_axis) {
                    return Err(limit());
                }
                root.refine_for_signs(400);
                let (lo, up) = root.isolator();
                let uu = I::new(lo.clone(), up.clone());
                let vv = enclose(&c0, &uu)
                    .div(&enclose(&c1, &uu))
                    .ok_or_else(limit)?
                    .neg();
                let x = |c: &R| I::exact(c.clone());
                let s1: [I; 3] = std::array::from_fn(|i| {
                    x(&o1[i]).add(&uu.mul(&x(&u[i]))).add(&vv.mul(&x(&v[i])))
                });
                let wv: [I; 3] = std::array::from_fn(|i| s1[i].sub(&x(&o2[i])));
                let h = wv
                    .iter()
                    .zip(a2)
                    .fold(I::exact_f64(0.0), |acc, (wi, ai)| acc.add(&wi.mul(&x(ai))))
                    .div(&x(&aa2))
                    .ok_or_else(limit)?;
                let m: [I; 3] = std::array::from_fn(|i| wv[i].sub(&h.mul(&x(&a2[i]))));
                let rho = m
                    .iter()
                    .fold(I::exact_f64(0.0), |acc, g| acc.add(&g.square()))
                    .sqrt();
                let qv = x(k)
                    .sub(
                        &wv.iter()
                            .fold(I::exact_f64(0.0), |acc, g| acc.add(&g.square())),
                    )
                    .sub(&x(&r22));
                let sigma = match qv.sign() {
                    Some(std::cmp::Ordering::Greater) => 1.0,
                    Some(std::cmp::Ordering::Less) => -1.0,
                    _ => return Err(limit()),
                };
                // s2 = o2 - sigma R2 M w / rho.
                let f = I::exact_f64(sigma)
                    .mul(&x(big2))
                    .div(&rho)
                    .ok_or_else(limit)?;
                let s2: [I; 3] = std::array::from_fn(|i| x(&o2[i]).sub(&f.mul(&m[i])));
                let gap: [I; 3] = std::array::from_fn(|i| s2[i].sub(&s1[i]));
                let norm = gap
                    .iter()
                    .fold(I::exact_f64(0.0), |acc, g| acc.add(&g.square()))
                    .sqrt();
                if norm.sign() != Some(std::cmp::Ordering::Greater) {
                    return Err(limit());
                }
                let away = k != &ks[0] && small2 > small1;
                let r = if away { x(small1).neg() } else { x(small1) };
                let point: [I; 3] = std::array::from_fn(|i| {
                    s1[i].add(&r.mul(&gap[i]).div(&norm).expect("a positive distance"))
                });
                found.push(Contact { spine: s1, point });
            }
        }
        if separated {
            return Ok(found);
        }
    }
    Err(limit())
}
