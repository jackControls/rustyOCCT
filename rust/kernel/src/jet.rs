//! Truncated Taylor series in one variable with enclosed coefficients
//! (S8d.2 of REVIEW_NOTES.md), and certified integration built on them.
//!
//! A `Jet` of order `n` holds `c[0..=n]`, enclosures of the Taylor
//! coefficients `f^(k)(t0) / k!` of a function about a base `t0`. With an
//! interval as the base, every coefficient encloses `f^(k)(xi) / k!` for
//! every `xi` in it (interval automatic differentiation), which bounds a
//! Taylor polynomial's Lagrange remainder. The operations are the standard
//! recurrences: sums and products term by term and by Cauchy products,
//! quotients and square roots by solving for the next coefficient, and the
//! transcendental functions through their derivatives' series integrated
//! term by term (`sin`/`cos` jointly, `atan` and `atan2` from `x y' - y x'
//! over x^2 + y^2`, `asin` from `f' / sqrt(1 - f^2)`), their constant terms
//! the certified values.
//!
//! `integrate` encloses `integral_a^b g(t) dt`: on each piece the Taylor
//! polynomial about its midpoint, integrated exactly, plus the remainder
//! `sup |c_{n+1}| 2 h^{n+2} / (n + 2)` from the piece's interval jet; a
//! piece whose remainder exceeds its share of the requested width is
//! bisected, down to a depth limit.
use crate::certified::Real;

/// A truncated Taylor series with enclosed coefficients.
#[derive(Debug, Clone)]
pub(crate) struct Jet<T> {
    pub(crate) c: Vec<T>,
}

fn zero<T: Real>() -> T {
    T::exact_f64(0.0)
}

impl<T: Real> Jet<T> {
    pub(crate) fn order(&self) -> usize {
        self.c.len() - 1
    }

    pub(crate) fn constant(x: T, n: usize) -> Self {
        let mut c = vec![zero::<T>(); n + 1];
        c[0] = x;
        Self { c }
    }

    /// The variable itself about a base: `t0 + s`.
    pub(crate) fn variable(t0: T, n: usize) -> Self {
        let mut c = vec![zero::<T>(); n + 1];
        c[0] = t0;
        if n >= 1 {
            c[1] = T::exact_f64(1.0);
        }
        Self { c }
    }

    pub(crate) fn add(&self, o: &Self) -> Self {
        Self {
            c: self.c.iter().zip(&o.c).map(|(a, b)| a.add(b)).collect(),
        }
    }

    pub(crate) fn sub(&self, o: &Self) -> Self {
        Self {
            c: self.c.iter().zip(&o.c).map(|(a, b)| a.sub(b)).collect(),
        }
    }

    pub(crate) fn neg(&self) -> Self {
        Self {
            c: self.c.iter().map(Real::neg).collect(),
        }
    }

    pub(crate) fn scale(&self, k: &T) -> Self {
        Self {
            c: self.c.iter().map(|a| a.mul(k)).collect(),
        }
    }

    pub(crate) fn add_constant(&self, k: &T) -> Self {
        let mut out = self.clone();
        out.c[0] = out.c[0].add(k);
        out
    }

    pub(crate) fn mul(&self, o: &Self) -> Self {
        let n = self.order();
        let c = (0..=n).map(|k| T::convolve(&self.c, &o.c, k)).collect();
        Self { c }
    }

    pub(crate) fn square(&self) -> Self {
        self.mul(self)
    }

    /// `self / o`: `q_k = (a_k - sum_{i<k} q_i b_{k-i}) / b_0`.
    pub(crate) fn div(&self, o: &Self) -> Option<Self> {
        let n = self.order();
        let mut q: Vec<T> = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let acc = if k == 0 {
                self.c[0].clone()
            } else {
                self.c[k].sub(&T::convolve(&q, &o.c[1..], k - 1))
            };
            q.push(acc.div(&o.c[0])?);
        }
        Some(Self { c: q })
    }

    /// `sqrt(self)`, `s_k = (a_k - sum_{0<i<k} s_i s_{k-i}) / (2 s_0)`.
    pub(crate) fn sqrt(&self) -> Option<Self> {
        let n = self.order();
        let s0 = self.c[0].sqrt();
        let two_s0 = s0.mul(&T::exact_f64(2.0));
        let mut s = vec![s0];
        for k in 1..=n {
            let acc = if k < 2 {
                self.c[k].clone()
            } else {
                self.c[k].sub(&T::convolve(&s[1..], &s[1..], k - 2))
            };
            s.push(acc.div(&two_s0)?);
        }
        Some(Self { c: s })
    }

    /// The series' derivative, one order shorter (its last coefficient
    /// zero).
    pub(crate) fn derivative(&self) -> Self {
        let n = self.order();
        let mut c: Vec<T> = (1..=n)
            .map(|k| self.c[k].mul(&T::exact_f64(k as f64)))
            .collect();
        c.push(zero::<T>());
        Self { c }
    }

    /// The antiderivative of `d` with constant term `c0`.
    fn integral(d: &Self, c0: T) -> Self {
        let n = d.order();
        let mut c = vec![c0];
        for k in 1..=n {
            c.push(
                d.c[k - 1]
                    .div(&T::exact_f64(k as f64))
                    .expect("a nonzero index"),
            );
        }
        Self { c }
    }

    /// `(cos self, sin self)`: `s' = c f'`, `c' = -s f'`, solved term by
    /// term from the certified values at the base.
    pub(crate) fn cos_sin(&self) -> (Self, Self) {
        let n = self.order();
        let (c0, s0) = T::cos_sin(&self.c[0]);
        let mut cs = vec![c0];
        let mut sn = vec![s0];
        // f' coefficients: (k + 1) f_{k+1}.
        let d: Vec<T> = (0..n)
            .map(|k| self.c[k + 1].mul(&T::exact_f64((k + 1) as f64)))
            .collect();
        for k in 1..=n {
            // k s_k = sum_{i<k} c_i d_{k-1-i}; k c_k = -sum s_i d_{k-1-i}.
            let (a, b) = (T::convolve(&cs, &d, k - 1), T::convolve(&sn, &d, k - 1));
            let kk = T::exact_f64(k as f64);
            sn.push(a.div(&kk).expect("a nonzero index"));
            cs.push(b.neg().div(&kk).expect("a nonzero index"));
        }
        (Self { c: cs }, Self { c: sn })
    }

    /// `e^self`: `E' = E f'`, solved term by term from the certified value
    /// at the base.
    pub(crate) fn exp(&self) -> Self {
        let n = self.order();
        let mut e = vec![T::exp(&self.c[0])];
        let d: Vec<T> = (0..n)
            .map(|k| self.c[k + 1].mul(&T::exact_f64((k + 1) as f64)))
            .collect();
        for k in 1..=n {
            let a = T::convolve(&e, &d, k - 1);
            e.push(a.div(&T::exact_f64(k as f64)).expect("a nonzero index"));
        }
        Self { c: e }
    }

    /// `(cosh self, sinh self)` from `e^self` and `e^-self`.
    pub(crate) fn cosh_sinh(&self) -> (Self, Self) {
        let (a, b) = (self.exp(), self.neg().exp());
        let half = T::exact_f64(0.5);
        (a.add(&b).scale(&half), a.sub(&b).scale(&half))
    }

    /// `atan2(y, x)` with the constant term `theta0` (a branch chosen by the
    /// caller, an enclosure of the angle at the base): its derivative is
    /// `(x y' - y x') / (x^2 + y^2)`.
    pub(crate) fn atan2(y: &Self, x: &Self, theta0: T) -> Option<Self> {
        let num = x.mul(&y.derivative()).sub(&y.mul(&x.derivative()));
        let den = x.square().add(&y.square());
        Some(Self::integral(&num.div(&den)?, theta0))
    }

    #[cfg(test)]
    /// `asin(self)` from `asin' = f' / sqrt(1 - f^2)`, its constant term
    /// `atan2(f0, sqrt(1 - f0^2))`.
    pub(crate) fn asin(&self) -> Option<Self> {
        let one = Self::constant(T::exact_f64(1.0), self.order());
        let root = one.sub(&self.square()).sqrt()?;
        let a0 = T::atan2(&self.c[0], &root.c[0])?;
        Some(Self::integral(&self.derivative().div(&root)?, a0))
    }
}

/// An enclosure of `integral_a^b g(t) dt` for `a < b`: `g` maps the jet of
/// the variable (about a piece's midpoint, or over the whole piece) to the
/// integrand's jet of the same order. Each piece integrates its Taylor
/// polynomial of order `order` about a midpoint `m` exactly over `[lo - m,
/// hi - m]` and bounds the remainder by its interval jet's next coefficient
/// times `(|lo - m|^{n+2} + |hi - m|^{n+2}) / (n + 2)`; pieces whose
/// remainder is wider than their share of `width` are bisected, at most
/// `depth` times. `None` when an evaluation fails (a
/// division by an enclosure of zero).
#[cfg(test)]
pub(crate) fn integrate<T: Real>(
    g: &dyn Fn(&Jet<T>) -> Option<Jet<T>>,
    a: f64,
    b: f64,
    order: usize,
    width: f64,
    depth: usize,
) -> Option<T> {
    let many = |t: &Jet<T>| g(t).map(|j| vec![j]);
    integrate_many(&many, 1, a, b, order, width, depth, false)?.pop()
}

/// Integrands of one variable's jet, evaluated together.
pub(crate) type Integrands<'a, T> = &'a dyn Fn(&Jet<T>) -> Option<Vec<Jet<T>>>;

/// [`integrate`] of `n` integrands at once from one evaluation per piece
/// (the jets they share computed once): a piece is bisected until every
/// integrand's remainder is within its share, of `width` itself or, when
/// `relative`, of `width` times the integrand's scale (its largest value at
/// eight points, at least one: mass moments, not signs). A bisected piece's
/// next coefficients, enclosed over it, hold over its halves too: a half
/// whose remainder they already bound within its share is integrated
/// without jets over it of its own (S9f.3's loops: most of the last
/// bisection's halves, a tenth of the integrals' time).
#[allow(clippy::too_many_arguments)]
pub(crate) fn integrate_many<T: Real>(
    g: Integrands<'_, T>,
    n: usize,
    a: f64,
    b: f64,
    order: usize,
    width: f64,
    depth: usize,
    relative: bool,
) -> Option<Vec<T>> {
    let mut totals = vec![zero::<T>(); n];
    // Each piece with its parent's next coefficients' magnitudes, if known.
    type Sizes = Option<std::rc::Rc<Vec<f64>>>;
    let mut stack: Vec<(f64, f64, usize, Sizes)> = vec![(a, b, 0usize, None)];
    let pow = |x: &T, k: usize| (0..k).fold(T::exact_f64(1.0), |p, _| p.mul(x));
    let mag = |x: &T| {
        let (xl, xh) = x.bounds_f64();
        T::exact_f64(xl.abs().max(xh.abs()))
    };
    // Each integrand's scale when the width asked for is relative to it: its
    // largest value at eight points, at least one.
    let mut scales = vec![1.0f64; n];
    if relative {
        for k in 0..8 {
            let t = a + (b - a) * (f64::from(k) + 0.5) / 8.0;
            if let Some(at) = g(&Jet::variable(T::exact_f64(t), 0)) {
                for (scale, jet) in scales.iter_mut().zip(&at) {
                    let (lo, hi) = jet.c[0].bounds_f64();
                    let size = lo.abs().max(hi.abs());
                    if size.is_finite() {
                        *scale = scale.max(size);
                    }
                }
            }
        }
    }
    while let Some((lo, hi, level, inherited)) = stack.pop() {
        let mid = 0.5 * lo + 0.5 * hi;
        let (l, m, h) = (T::exact_f64(lo), T::exact_f64(mid), T::exact_f64(hi));
        let (s0, s1) = (l.sub(&m), h.sub(&m));
        let reach = pow(&mag(&s0), order + 2).add(&pow(&mag(&s1), order + 2));
        let share = width * (hi - lo) / (b - a);
        // A coefficient's magnitude times the piece's reach over `order + 2`,
        // rounded up.
        let bound_of = |size: f64| {
            T::exact_f64(size)
                .mul(&reach)
                .div(&T::exact_f64((order + 2) as f64))
                .map_or(f64::INFINITY, |x| x.bounds_f64().1)
        };
        // A NaN bound never settles.
        let settled = |bounds: &[f64]| {
            bounds
                .iter()
                .zip(&scales)
                .all(|(bound, scale)| *bound <= share * scale)
        };
        let from_parent: Option<Vec<f64>> = inherited
            .as_ref()
            .map(|sizes| sizes.iter().map(|x| bound_of(*x)).collect())
            .filter(|bounds: &Vec<f64>| settled(bounds));
        let (bounds, sizes) = match from_parent {
            Some(bounds) => (bounds, None),
            None => {
                // The remainders: the next coefficients over the piece
                // (unbounded where the integrands' jets are not defined
                // over all of it).
                match g(&Jet::variable(l.union(&h), order + 1)) {
                    Some(over) if over.len() == n => {
                        let sizes: Vec<f64> = over
                            .iter()
                            .map(|j| {
                                let (nlo, nhi) = j.c[order + 1].bounds_f64();
                                nlo.abs().max(nhi.abs())
                            })
                            .collect();
                        (sizes.iter().map(|x| bound_of(*x)).collect(), Some(sizes))
                    }
                    _ => (vec![f64::INFINITY; n], None),
                }
            }
        };
        if !settled(&bounds) && level < depth && lo < mid && mid < hi {
            let sizes = sizes
                .filter(|x| x.iter().all(|y| y.is_finite()))
                .map(std::rc::Rc::new);
            stack.push((lo, mid, level + 1, sizes.clone()));
            stack.push((mid, hi, level + 1, sizes));
            continue;
        }
        if !bounds.iter().all(|bound| bound.is_finite()) {
            return None;
        }
        let at = g(&Jet::variable(m, order))?;
        if at.len() != n {
            return None;
        }
        for ((total, jet), bound) in totals.iter_mut().zip(&at).zip(&bounds) {
            let mut piece = zero::<T>();
            for k in 0..=order {
                let span = pow(&s1, k + 1).sub(&pow(&s0, k + 1));
                piece = piece.add(&jet.c[k].mul(&span).div(&T::exact_f64((k + 1) as f64))?);
            }
            let remainder = T::exact_f64(-bound).union(&T::exact_f64(*bound));
            *total = total.add(&piece).add(&remainder);
        }
    }
    Some(totals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::certified::Fast;

    fn contains(x: &Fast, v: f64) -> bool {
        let (lo, hi) = x.bounds_f64();
        lo <= v && v <= hi
    }

    /// Jets of known functions agree with their Taylor coefficients.
    #[test]
    fn jets_match_known_series() {
        let t = Jet::variable(Fast::exact_f64(0.3), 6);
        let (c, s) = t.cos_sin();
        // d^k/dt^k sin at 0.3 over k!.
        let fact = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0];
        for k in 0..=6 {
            let ds = [0.3f64.sin(), 0.3f64.cos(), -0.3f64.sin(), -0.3f64.cos()][k % 4] / fact[k];
            let dc = [0.3f64.cos(), -0.3f64.sin(), -0.3f64.cos(), 0.3f64.sin()][k % 4] / fact[k];
            assert!((s.c[k].bounds_f64().0 - ds).abs() < 1e-14, "sin {k}");
            assert!((c.c[k].bounds_f64().0 - dc).abs() < 1e-14, "cos {k}");
        }
        // atan2(sin t, cos t) = t.
        let a = Jet::atan2(&s, &c, Fast::atan2(&s.c[0], &c.c[0]).unwrap()).unwrap();
        assert!(contains(&a.c[0], 0.3) && contains(&a.c[1], 1.0));
        for k in 2..=6 {
            assert!(a.c[k].bounds_f64().1.abs() < 1e-12, "{k}");
        }
        // asin(sin t) = t.
        let b = s.asin().unwrap();
        assert!(contains(&b.c[0], 0.3) && contains(&b.c[1], 1.0));
        // sqrt(t^2) = t for t > 0; 1 / (1 + t) = sum (-t)^k about 0.
        let q = t.square().sqrt().unwrap();
        assert!(contains(&q.c[0], 0.3) && contains(&q.c[1], 1.0));
        let one = Jet::constant(Fast::exact_f64(1.0), 6);
        let u = Jet::variable(Fast::exact_f64(0.0), 6);
        let inv = one.div(&one.add(&u)).unwrap();
        for k in 0..=6 {
            assert!(
                contains(&inv.c[k], if k % 2 == 0 { 1.0 } else { -1.0 }),
                "{k}"
            );
        }
    }

    /// `exp` and its jet agree with the series; cosh^2 - sinh^2 = 1.
    #[test]
    fn exponentials_are_enclosed() {
        for x in [-3.0, -0.4, 0.0, 0.7, 2.5, 20.0] {
            let e = Fast::exp(&Fast::exact_f64(x));
            let (lo, hi) = e.bounds_f64();
            assert!(
                lo <= x.exp() * (1.0 + 1e-15) && x.exp() * (1.0 - 1e-15) <= hi,
                "{x}"
            );
            assert!(hi - lo <= 1e-12 * x.exp().max(1.0), "{x}: {lo} {hi}");
        }
        let t = Jet::variable(Fast::exact_f64(0.6), 5);
        let (c, s) = t.cosh_sinh();
        let one = c.square().sub(&s.square());
        assert!(contains(&one.c[0], 1.0));
        for k in 1..=5 {
            assert!(one.c[k].bounds_f64().1.abs() < 1e-12, "{k}");
        }
    }

    /// Certified integrals contain the exact values and are narrow.
    #[test]
    fn integrals_enclose_the_exact_values() {
        // integral_0^pi sin t dt = 2.
        let g = |t: &Jet<Fast>| Some(t.cos_sin().1);
        let v = integrate(&g, 0.0, std::f64::consts::PI, 6, 1e-12, 30).unwrap();
        let (lo, hi) = v.bounds_f64();
        assert!(lo <= 2.0 && 2.0 <= hi && hi - lo < 1e-10, "{lo} {hi}");
        // integral_0^1 1 / (1 + t^2) dt = pi / 4.
        let g = |t: &Jet<Fast>| {
            let one = Jet::constant(Fast::exact_f64(1.0), t.order());
            one.div(&one.add(&t.square()))
        };
        let v = integrate(&g, 0.0, 1.0, 6, 1e-12, 30).unwrap();
        let (lo, hi) = v.bounds_f64();
        let q = std::f64::consts::FRAC_PI_4;
        assert!(lo <= q && q <= hi && hi - lo < 1e-10, "{lo} {hi}");
    }
}
