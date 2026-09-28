//! Tensor Bernstein polynomials on `[0, 1]²` in the `Fast` tier, for the
//! rational cells' derivative numerators and the normals' cones (T-b). A
//! polynomial's values lie in the convex hull of its coefficients, so the
//! largest coefficient bounds it; products and derivatives are computed on
//! the coefficients, which keeps exact cancellations (a rational extrusion's
//! `S_uv`) cancelled up to rounding.
use crate::certified::{Fast, Real};

fn c(x: f64) -> Fast {
    Fast::exact_f64(x)
}

/// `C(n, k)`, exact in binary64 for `n <= 56`.
fn binomial(n: usize, k: usize) -> f64 {
    let mut x = 1.0f64;
    for i in 0..k {
        x = x * (n - i) as f64 / (i + 1) as f64;
    }
    x.round()
}

/// Coefficients of degrees `d` in `(u, v)`, `c[i (d[1] + 1) + j]`.
#[derive(Debug, Clone)]
pub(super) struct Poly {
    d: [usize; 2],
    c: Vec<Fast>,
}

impl Poly {
    pub(super) fn new(d: [usize; 2], c: Vec<Fast>) -> Self {
        debug_assert_eq!(c.len(), (d[0] + 1) * (d[1] + 1));
        Self { d, c }
    }

    fn zero(d: [usize; 2]) -> Self {
        Self::new(d, vec![c(0.0); (d[0] + 1) * (d[1] + 1)])
    }

    fn at(&self, i: usize, j: usize) -> &Fast {
        &self.c[i * (self.d[1] + 1) + j]
    }

    pub(super) fn coefficients(&self) -> &[Fast] {
        &self.c
    }

    /// The value at an enclosed `(u, v)` in `[0, 1]²`, by de Casteljau.
    pub(super) fn value(&self, u: &Fast, v: &Fast) -> Fast {
        let casteljau = |row: Vec<Fast>, t: &Fast| {
            let s = c(1.0).sub(t);
            let mut row = row;
            while row.len() > 1 {
                row = row
                    .windows(2)
                    .map(|w| w[0].mul(&s).add(&w[1].mul(t)))
                    .collect();
            }
            row[0]
        };
        let rows: Vec<Fast> = (0..=self.d[0])
            .map(|i| casteljau((0..=self.d[1]).map(|j| *self.at(i, j)).collect(), v))
            .collect();
        casteljau(rows, u)
    }

    /// The same polynomial of higher degrees.
    fn elevated(&self, d: [usize; 2]) -> Self {
        let mut out = self.clone();
        while out.d[0] < d[0] {
            let n = out.d[0];
            let mut next = Self::zero([n + 1, out.d[1]]);
            for i in 0..=n + 1 {
                let f = c(i as f64)
                    .div(&c((n + 1) as f64))
                    .expect("a positive degree");
                let g = c(1.0).sub(&f);
                for j in 0..=out.d[1] {
                    let left = if i > 0 {
                        f.mul(out.at(i - 1, j))
                    } else {
                        c(0.0)
                    };
                    let right = if i <= n { g.mul(out.at(i, j)) } else { c(0.0) };
                    next.c[i * (out.d[1] + 1) + j] = left.add(&right);
                }
            }
            out = next;
        }
        while out.d[1] < d[1] {
            let n = out.d[1];
            let mut next = Self::zero([out.d[0], n + 1]);
            for j in 0..=n + 1 {
                let f = c(j as f64)
                    .div(&c((n + 1) as f64))
                    .expect("a positive degree");
                let g = c(1.0).sub(&f);
                for i in 0..=out.d[0] {
                    let left = if j > 0 {
                        f.mul(out.at(i, j - 1))
                    } else {
                        c(0.0)
                    };
                    let right = if j <= n { g.mul(out.at(i, j)) } else { c(0.0) };
                    next.c[i * (n + 2) + j] = left.add(&right);
                }
            }
            out = next;
        }
        out
    }

    fn combine(&self, o: &Self, f: impl Fn(&Fast, &Fast) -> Fast) -> Self {
        let d = [self.d[0].max(o.d[0]), self.d[1].max(o.d[1])];
        let (a, b) = (self.elevated(d), o.elevated(d));
        Self::new(d, a.c.iter().zip(&b.c).map(|(x, y)| f(x, y)).collect())
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        self.combine(o, |x, y| x.add(y))
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        self.combine(o, |x, y| x.sub(y))
    }

    pub(super) fn scale(&self, k: f64) -> Self {
        Self::new(self.d, self.c.iter().map(|x| x.mul(&c(k))).collect())
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        let d = [self.d[0] + o.d[0], self.d[1] + o.d[1]];
        let scaled = |p: &Self| -> Vec<Fast> {
            let mut out = Vec::with_capacity(p.c.len());
            for i in 0..=p.d[0] {
                for j in 0..=p.d[1] {
                    let k = binomial(p.d[0], i) * binomial(p.d[1], j);
                    out.push(p.at(i, j).mul(&c(k)));
                }
            }
            out
        };
        let (a, b) = (scaled(self), scaled(o));
        let mut out = Self::zero(d);
        for i in 0..=self.d[0] {
            for j in 0..=self.d[1] {
                let x = &a[i * (self.d[1] + 1) + j];
                for k in 0..=o.d[0] {
                    for l in 0..=o.d[1] {
                        let at = (i + k) * (d[1] + 1) + j + l;
                        out.c[at] = out.c[at].add(&x.mul(&b[k * (o.d[1] + 1) + l]));
                    }
                }
            }
        }
        for i in 0..=d[0] {
            for j in 0..=d[1] {
                let k = binomial(d[0], i) * binomial(d[1], j);
                let at = i * (d[1] + 1) + j;
                out.c[at] = out.c[at].div(&c(k)).expect("a positive binomial");
            }
        }
        out
    }

    /// The derivative in `u` (axis 0) or `v` (axis 1) on `[0, 1]`; a constant
    /// direction's is zero.
    pub(super) fn derivative(&self, axis: usize) -> Self {
        let n = self.d[axis];
        if n == 0 {
            return Self::zero(self.d);
        }
        let mut d = self.d;
        d[axis] = n - 1;
        let mut out = Self::zero(d);
        let k = c(n as f64);
        for i in 0..=d[0] {
            for j in 0..=d[1] {
                let (a, b) = if axis == 0 {
                    (self.at(i + 1, j), self.at(i, j))
                } else {
                    (self.at(i, j + 1), self.at(i, j))
                };
                out.c[i * (d[1] + 1) + j] = a.sub(b).mul(&k);
            }
        }
        out
    }
}

/// A vector of three polynomials.
pub(super) type Poly3 = [Poly; 3];

pub(super) fn largest3(p: &Poly3) -> f64 {
    // The coefficient vectors' lengths bound the vector polynomial.
    let (a, b, cc) = (&p[0].c, &p[1].c, &p[2].c);
    let n = a.len().min(b.len()).min(cc.len());
    if a.len() != b.len() || b.len() != cc.len() {
        return f64::INFINITY;
    }
    (0..n)
        .map(|k| {
            a[k].square()
                .add(&b[k].square())
                .add(&cc[k].square())
                .sqrt()
                .bounds_f64()
                .1
        })
        .fold(0.0, f64::max)
}

pub(super) fn mul3(p: &Poly3, s: &Poly) -> Poly3 {
    [p[0].mul(s), p[1].mul(s), p[2].mul(s)]
}

pub(super) fn sub33(a: &Poly3, b: &Poly3) -> Poly3 {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

pub(super) fn add33(a: &Poly3, b: &Poly3) -> Poly3 {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}

pub(super) fn scale3(a: &Poly3, k: f64) -> Poly3 {
    [a[0].scale(k), a[1].scale(k), a[2].scale(k)]
}

pub(super) fn derivative3(a: &Poly3, axis: usize) -> Poly3 {
    [
        a[0].derivative(axis),
        a[1].derivative(axis),
        a[2].derivative(axis),
    ]
}

pub(super) fn cross3(a: &Poly3, b: &Poly3) -> Poly3 {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}
