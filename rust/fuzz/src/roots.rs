//! Structured complete-root oracle from products of known real/complex factors.
//! Query signs reduce directly in Q(sqrt(d)); no Sturm sequence is used here.
use super::{byte, zero};
use crate::splines::rat;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::polynomial::{Polynomial, RealRoots};
use rusty_occt::{Error, ScalarInterval};
use std::cmp::Ordering::{self, Equal, Greater, Less};

#[derive(Clone)]
struct Root {
    rational: Option<i32>,
    square: i32,
    positive: bool,
    multiplicity: usize,
}
impl Root {
    fn compare(&self, x: &R) -> Ordering {
        if let Some(n) = self.rational {
            return R::from_integer(BigInt::from(n)).cmp(x);
        }
        if self.positive {
            if x < &zero() {
                Greater
            } else {
                R::from_integer(BigInt::from(self.square)).cmp(&(x * x))
            }
        } else if x > &zero() {
            Less
        } else {
            (x * x).cmp(&R::from_integer(BigInt::from(self.square)))
        }
    }
    fn sign_at(&self, q: &[f64]) -> Ordering {
        if let Some(n) = self.rational {
            let r = R::from_integer(BigInt::from(n));
            return q
                .iter()
                .rev()
                .fold(zero(), |a, &c| a * &r + rat(c))
                .cmp(&zero());
        }
        let (mut a, mut b) = (zero(), zero());
        let d = R::from_integer(BigInt::from(self.square));
        for &c in q.iter().rev() {
            (a, b) = (b * &d + rat(c), a);
        }
        if b == zero() {
            a.cmp(&zero())
        } else {
            let result = self.compare(&(-a / &b));
            if b < zero() {
                result.reverse()
            } else {
                result
            }
        }
    }
}

pub(crate) fn interval(bounds: ScalarInterval, compare: impl Fn(&R) -> Ordering) {
    let (lo, hi) = (bounds.lower(), bounds.upper());
    if lo == hi {
        assert_eq!(compare(&rat(lo)), Equal);
    } else {
        assert_eq!(lo.next_up(), hi);
        assert_eq!(compare(&rat(lo)), Greater);
        assert_eq!(compare(&rat(hi)), Less);
    }
}

pub fn check_roots(data: &[u8]) {
    let n = 1 + usize::from(byte(data, 0) % 25);
    let mut coefficients = vec![1i64];
    // Exact ascending order, independent of any returned root bounds.
    let mut expected: Vec<Root> = [
        (None, 5, false),
        (Some(-2), 0, false),
        (None, 3, false),
        (None, 2, false),
        (Some(-1), 0, false),
        (Some(0), 0, false),
        (Some(1), 0, false),
        (None, 2, true),
        (None, 3, true),
        (Some(2), 0, true),
        (None, 5, true),
    ]
    .into_iter()
    .map(|(rational, square, positive)| Root {
        rational,
        square,
        positive,
        multiplicity: 0,
    })
    .collect();
    let mut degree = 0;
    while degree < n {
        let v = byte(data, 2 + degree);
        let factor = if degree + 2 <= n && v % 3 != 0 {
            let d = if v % 3 == 1 {
                [2, 3, 5][usize::from(v / 3) % 3]
            } else {
                -1
            };
            if d > 0 {
                for root in &mut expected {
                    if root.rational.is_none() && root.square == d {
                        root.multiplicity += 1;
                    }
                }
            }
            vec![-i64::from(d), 0, 1]
        } else {
            let root = i32::from(v % 5) - 2;
            expected
                .iter_mut()
                .find(|r| r.rational == Some(root))
                .unwrap()
                .multiplicity += 1;
            vec![-i64::from(root), 1]
        };
        degree += factor.len() - 1;
        let mut next = vec![0; coefficients.len() + factor.len() - 1];
        for (i, a) in coefficients.iter().enumerate() {
            for (j, b) in factor.iter().enumerate() {
                next[i + j] += a * b;
            }
        }
        coefficients = next;
    }
    let scale = 2f64.powi((i32::from(byte(data, 1)) - 128) * 5)
        * if byte(data, 27) & 1 == 0 { 1. } else { -1. };
    let p: Vec<f64> = coefficients.iter().map(|&x| x as f64 * scale).collect();
    for (&represented, &integer) in p.iter().zip(&coefficients) {
        assert_eq!(
            rat(represented),
            R::from_integer(BigInt::from(integer)) * rat(scale)
        );
    }
    let polynomial = Polynomial::new(&p).unwrap();
    let clipped = byte(data, 28) & 1 == 1;
    expected.retain(|r| {
        r.multiplicity > 0
            && (!clipped || (r.compare(&rat(-2.)) != Less && r.compare(&rat(2.)) != Greater))
    });
    let RealRoots::Finite(actual) = (if clipped {
        polynomial.roots_in(-2., 2.)
    } else {
        polynomial.real_roots()
    })
    .unwrap() else {
        panic!("nonzero product returned all roots");
    };
    assert_eq!(actual.len(), expected.len());
    let q: Vec<f64> = (0..=usize::from(byte(data, 29) % 25))
        .map(|i| {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 30 + 8 * i + j)
            })))
        })
        .collect();
    let query = Polynomial::new(&q);
    let finite = q.iter().all(|x| x.is_finite());
    if !finite {
        assert!(matches!(query, Err(Error::NonFinite(_))));
    }
    for (index, (root, expected)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(root.multiplicity(), expected.multiplicity);
        interval(root.bounds().unwrap(), |x| expected.compare(x));
        assert_eq!(root.sign_at(&polynomial), Equal);
        // These one- or two-root equations independently describe the known
        // factor. Equality must survive different defining polynomials and
        // multiplicities; ordering uses the independently ordered factor list.
        let factor = if let Some(n) = expected.rational {
            vec![-f64::from(n), 1.]
        } else {
            vec![-f64::from(expected.square), 0., 1.]
        };
        let RealRoots::Finite(factor_roots) =
            Polynomial::new(&factor).unwrap().real_roots().unwrap()
        else {
            panic!("known factor is not zero");
        };
        let known = &factor_roots[usize::from(expected.rational.is_none() && expected.positive)];
        assert_eq!(root.compare_root(known), Equal);
        assert_eq!(known.compare_root(root), Equal);
        let other = usize::from(byte(data, 238 + index)) % actual.len();
        assert_eq!(known.compare_root(&actual[other]), index.cmp(&other));
        assert_eq!(actual[other].compare_root(known), other.cmp(&index));
        if finite {
            assert_eq!(root.sign_at(query.as_ref().unwrap()), expected.sign_at(&q));
        }
    }
    if let Some(&probe) = q.first() {
        if finite {
            assert_eq!(
                polynomial.sign(probe).unwrap(),
                p.iter()
                    .rev()
                    .fold(zero(), |a, &c| a * rat(probe) + rat(c))
                    .cmp(&zero())
            );
        } else if !probe.is_finite() {
            assert!(matches!(polynomial.sign(probe), Err(Error::NonFinite(_))));
        }
    }
    assert_eq!(Polynomial::new(&[0.]).unwrap().sign(0.).unwrap(), Equal);
    assert!(matches!(
        Polynomial::new(&vec![1.; 27]),
        Err(Error::LimitExceeded(_))
    ));
}
