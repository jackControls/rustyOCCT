//! Deterministic representation changes shared by tests and the native bridge.
use num_rational::BigRational as R;
use rusty_occt::ExactBSplineCurve3;

pub fn refined(curve: &ExactBSplineCurve3) -> (ExactBSplineCurve3, [R; 2]) {
    let a = &curve.domain()[0];
    let b = curve.knots().iter().find(|k| *k > a).unwrap();
    let u = a + (b - a) / R::from_integer(3.into());
    let v = a + (b - a) * R::new(2.into(), 3.into());
    let fine = curve
        .refined(&[(v.clone(), 1), (u.clone(), curve.degree()), (v.clone(), 0)])
        .unwrap();
    (fine, [u, v])
}

pub fn restored(curve: &ExactBSplineCurve3, cuts: &[R; 2]) -> ExactBSplineCurve3 {
    curve
        .remove_knot(&cuts[1], 0)
        .unwrap()
        .unwrap()
        .remove_knot(&cuts[0], 0)
        .unwrap()
        .unwrap()
}
