use rusty_occt::polynomial::{Polynomial, RealRoots, RootIsolationOptions};
use rusty_occt::Error;
use std::cmp::Ordering::{Equal, Greater, Less};

#[test]
fn roots_and_algebraic_queries_match_independent_continued_fraction_oracle() {
    fn number(word: &str) -> f64 {
        f64::from_bits(u64::from_str_radix(word, 16).unwrap())
    }
    let mut count = 0;
    for row in include_str!("../../fixtures/real-roots.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let mut w = row.split_whitespace();
        let name = w.next().unwrap();
        let np = w.next().unwrap().parse::<usize>().unwrap();
        let nq = w.next().unwrap().parse::<usize>().unwrap();
        let lower = w.next().unwrap();
        let upper = w.next().unwrap();
        let p = Polynomial::new(
            &(0..np)
                .map(|_| number(w.next().unwrap()))
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let q = Polynomial::new(
            &(0..nq)
                .map(|_| number(w.next().unwrap()))
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let found = if lower == "*" {
            p.real_roots()
        } else {
            p.roots_in(number(lower), number(upper))
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        let status = w.next().unwrap();
        match found {
            RealRoots::All => assert_eq!(status, "A", "{name}"),
            RealRoots::Finite(roots) => {
                assert_eq!(status, "N");
                assert_eq!(
                    roots.len(),
                    w.next().unwrap().parse::<usize>().unwrap(),
                    "{name}"
                );
                for r in &roots {
                    assert_eq!(
                        r.multiplicity(),
                        w.next().unwrap().parse::<usize>().unwrap(),
                        "{name}"
                    );
                    let status = w.next().unwrap();
                    if status == "U" {
                        assert!(
                            matches!(r.bounds(), Err(Error::Unrepresentable(_))),
                            "{name}"
                        );
                    } else {
                        assert_eq!(status, "E");
                        let lo = number(w.next().unwrap());
                        let hi = number(w.next().unwrap());
                        let bounds = r.bounds().unwrap();
                        assert_eq!((bounds.lower(), bounds.upper()), (lo, hi), "{name}");
                        assert_eq!(
                            r.compare(lo).unwrap(),
                            if lo == hi { Equal } else { Greater },
                            "{name}"
                        );
                        assert_eq!(
                            r.compare(hi).unwrap(),
                            if lo == hi { Equal } else { Less },
                            "{name}"
                        );
                    }
                    let sign = w.next().unwrap().parse::<i32>().unwrap().cmp(&0);
                    assert_eq!(r.sign_at(&q), sign, "{name}");
                    assert_eq!(r.sign_at(&p), Equal, "{name}");
                }
                if name.starts_with("cluster") {
                    assert_eq!(roots.len(), 2);
                    assert_eq!(roots[0].bounds().unwrap(), roots[1].bounds().unwrap());
                    assert_eq!(roots[0].sign_at(&q), Less);
                    assert_eq!(roots[1].sign_at(&q), Greater);
                }
            }
        }
        assert!(w.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 94);
}

#[test]
fn repeated_roots_and_signs_at_algebraic_arguments_are_exact() {
    // (x²-2)² (x-1)³
    let polynomial = Polynomial::new(&[-4., 12., -8., -8., 11., -1., -3., 1.]).unwrap();
    let RealRoots::Finite(roots) = polynomial.real_roots().unwrap() else {
        panic!()
    };
    assert_eq!(roots.len(), 3);
    assert_eq!(
        roots.iter().map(|r| r.multiplicity()).collect::<Vec<_>>(),
        vec![2, 3, 2]
    );
    for r in &roots {
        assert_eq!(r.sign_at(&polynomial), Equal);
    }
    let sqrt = Polynomial::new(&[-2., 0., 1.]).unwrap();
    assert_eq!(roots[0].sign_at(&sqrt), Equal);
    assert_eq!(roots[2].sign_at(&sqrt), Equal);
    assert_eq!(roots[1].sign_at(&sqrt), Less);
    assert_eq!(roots[1].compare(1.).unwrap(), Equal);
    assert_eq!(roots[0].compare(-1.5).unwrap(), Greater);
    assert_eq!(roots[0].compare(-1.4).unwrap(), Less);
    assert_eq!(roots[2].compare(1.4).unwrap(), Greater);
    assert_eq!(roots[2].compare(1.5).unwrap(), Less);
}

#[test]
fn closed_intervals_include_endpoints_once_and_resource_failure_is_atomic() {
    let p = Polynomial::new(&[0., -1., 0., 1.]).unwrap();
    for (a, b, expected) in [
        (-1., 1., vec![-1., 0., 1.]),
        (0., 0., vec![0.]),
        (0., 1., vec![0., 1.]),
        (0.25, 0.75, vec![]),
    ] {
        let RealRoots::Finite(roots) = p.roots_in(a, b).unwrap() else {
            panic!()
        };
        assert_eq!(roots.len(), expected.len());
        for (r, e) in roots.iter().zip(expected) {
            assert_eq!(r.compare(e).unwrap(), Equal);
        }
    }
    assert!(matches!(
        p.real_roots_with_options(RootIsolationOptions {
            max_subdivisions: 0
        }),
        Err(Error::ComputationLimit(_))
    ));
    assert!(matches!(p.roots_in(2., 1.), Err(Error::OutOfDomain(_))));
    assert!(Polynomial::new(&vec![0.; 27]).is_err());
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(Polynomial::new(&[x]), Err(Error::NonFinite(_))));
        assert!(matches!(p.roots_in(x, 1.), Err(Error::NonFinite(_))));
        assert!(matches!(p.sign(x), Err(Error::NonFinite(_))));
    }
    assert!(matches!(
        Polynomial::new(&[]).unwrap().real_roots().unwrap(),
        RealRoots::All
    ));
    let RealRoots::Finite(roots) = Polynomial::new(&[1.]).unwrap().real_roots().unwrap() else {
        panic!()
    };
    assert!(roots.is_empty());
}

#[test]
fn finite_coefficients_can_have_unrepresentable_roots() {
    let p = Polynomial::new(&[-f64::MAX, f64::from_bits(1)]).unwrap();
    let RealRoots::Finite(roots) = p.real_roots().unwrap() else {
        panic!()
    };
    assert_eq!(roots.len(), 1);
    assert_eq!(roots[0].compare(f64::MAX).unwrap(), Greater);
    assert!(matches!(roots[0].bounds(), Err(Error::Unrepresentable(_))));
    assert_eq!(roots[0].sign_at(&p), Equal);
}
