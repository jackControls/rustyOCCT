//! S10's budgets (`BUDGETS.md`): an operation under a budget of work units,
//! a deadline or a cancellation token either succeeds with the result it
//! has without one or is refused with `BudgetExhausted`, never a partial
//! body; work units are deterministic; another thread cancels.
use rusty_occt::budget::{Cause, Exhaustion};
use rusty_occt::identity::OperationId;
use rusty_occt::{Budget, Cancel, Context, Error, Frame3, Point2, Point3, Solid, Tolerance, Vec3};
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    Frame3::new(
        Point3::new(o[0], o[1], o[2]),
        Vec3::new(n[0], n[1], n[2]),
        Vec3::new(x[0], x[1], x[2]),
        tol(),
    )
    .unwrap()
}

/// `torus_hole_rod_booleans.rs`'s tilted frame: a whole torus of radii 3
/// and 1 and a rod of radius `2 + g` about its axis within rounding.
fn torus_and_rod(g: f64) -> (Solid, Solid) {
    let f = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let pi = std::f64::consts::PI;
    let torus = Solid::torus_with(OperationId(1), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
        .unwrap()
        .0;
    let at = |l: [f64; 3]| f.point(Point2::new(l[0], l[1]), l[2]).to_array();
    let (x, n) = (f.x().to_array(), f.normal().to_array());
    let rf = frame(at([0.0, 0.0, -1.5]), n, x);
    let rod = Solid::cylinder_with(OperationId(2), rf, 2.0 + g, 0.0, 3.0, tol())
        .unwrap()
        .0;
    (torus, rod)
}

/// Two boxes overlapping at a corner (the polyhedral engine).
fn boxes() -> (Solid, Solid) {
    let a = Solid::box_at_with(
        OperationId(1),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(2.0, 2.0, 2.0),
        tol(),
    )
    .unwrap()
    .0;
    let b = Solid::box_at_with(
        OperationId(2),
        Point3::new(1.0, 1.0, 1.0),
        Vec3::new(2.0, 2.0, 2.0),
        tol(),
    )
    .unwrap()
    .0;
    (a, b)
}

fn fuse_under(
    a: &Solid,
    b: &Solid,
    budget: &Budget,
) -> rusty_occt::Result<(Vec<Solid>, rusty_occt::history::History)> {
    a.fuse_in(&Context::new(OperationId(3)).with_budget(budget), b)
}

/// The refusal of a budget of `n` units.
fn refused_at(n: u64) -> Error {
    Error::BudgetExhausted(Exhaustion {
        cause: Cause::Work,
        stage: "",
        work: n + 1,
    })
}

fn same_refusal(got: &Error, want: &Error) -> bool {
    match (got, want) {
        (Error::BudgetExhausted(g), Error::BudgetExhausted(w)) => {
            g.cause == w.cause && g.work == w.work
        }
        _ => false,
    }
}

/// Without a budget, and under one that suffices, the same result; with
/// one unit fewer than it counts, refused at exactly one unit past the
/// budget; every count the same on a second run (a kept arrangement
/// charged again) and pinned.
#[test]
fn a_budget_that_suffices_keeps_the_result_and_one_unit_less_refuses() {
    // The counts are pinned: the same on every host (`BUDGETS.md`).
    for (name, (a, b), pinned) in [
        ("boxes", boxes(), 477),
        ("torus and rod", torus_and_rod(1e-6), 51_200),
    ] {
        let plain = a.fuse(OperationId(3), &b);
        let full = Budget::unlimited();
        let counted = fuse_under(&a, &b, &full);
        let w = full.spent();
        assert!(w > 0, "{name}: no work counted");
        assert_eq!(w, pinned, "{name}: the work units");
        assert_eq!(counted, plain, "{name}: an unlimited budget's result");
        let again = Budget::unlimited();
        assert_eq!(fuse_under(&a, &b, &again), plain, "{name}: again");
        assert_eq!(again.spent(), w, "{name}: the work again");
        let exact = Budget::work(w);
        assert_eq!(fuse_under(&a, &b, &exact), plain, "{name}: within {w}");
        assert_eq!(exact.spent(), w, "{name}");
        for n in [0, w / 2, w - 1] {
            let short = Budget::work(n);
            let got = fuse_under(&a, &b, &short).unwrap_err();
            assert!(
                same_refusal(&got, &refused_at(n)),
                "{name}: {got:?} under {n} of {w}"
            );
            assert_eq!(short.spent(), n + 1, "{name}: {n}");
            assert!(short.exhausted().unwrap().deterministic());
        }
    }
}

/// An exhausted budget refuses every later operation at once; clones
/// share one count.
#[test]
fn an_exhausted_budget_refuses_later_operations() {
    let (a, b) = boxes();
    let budget = Budget::work(1);
    let first = fuse_under(&a, &b, &budget).unwrap_err();
    let shared = budget.clone();
    let second = a
        .cut_in(&Context::new(OperationId(4)).with_budget(&shared), &b)
        .unwrap_err();
    assert_eq!(first, second);
    assert_eq!(shared.spent(), 2);
    // The validation of a body, under `Budget::run`.
    let topology = a.topology();
    assert_eq!(budget.run(|| topology.validate(tol())), Err(first));
}

/// A deadline already passed refuses for wall time, distinctly.
#[test]
fn a_passed_deadline_refuses_for_wall_time() {
    let (a, b) = boxes();
    let budget = Budget::unlimited().with_deadline(Instant::now());
    match fuse_under(&a, &b, &budget) {
        Err(Error::BudgetExhausted(e)) => {
            assert_eq!(e.cause, Cause::WallTime);
            assert!(!e.deterministic());
        }
        other => panic!("{other:?}"),
    }
}

/// Another thread cancels a Boolean that takes about a second in release
/// (the tilted torus and a rod crossing its inner equator by `1e-6`): the
/// operation returns `Cancelled` promptly.
#[test]
fn another_thread_cancels_a_long_boolean_promptly() {
    let (torus, rod) = torus_and_rod(1e-6);
    let cancel = Cancel::new();
    let budget = Budget::unlimited().with_cancel(&cancel);
    let canceller = {
        let cancel = cancel.clone();
        std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(50));
            cancel.cancel();
            Instant::now()
        })
    };
    let start = Instant::now();
    let out = std::thread::scope(|s| {
        s.spawn(|| fuse_under(&torus, &rod, &budget))
            .join()
            .expect("no panic")
    });
    let returned = Instant::now();
    let cancelled = canceller.join().unwrap();
    match out {
        Err(Error::BudgetExhausted(e)) => assert_eq!(e.cause, Cause::Cancelled, "{e:?}"),
        other => panic!(
            "not cancelled after {:?}: {:?}",
            start.elapsed(),
            other.map(|x| x.0.len())
        ),
    }
    let latency = returned.saturating_duration_since(cancelled);
    eprintln!(
        "cancelled after {:?}, returned {latency:?} later",
        cancelled - start
    );
    assert!(
        latency < Duration::from_millis(500),
        "cancelled after {latency:?}"
    );
    // The arrangement refused is not kept: the same operation unlimited
    // evaluates afterwards on this thread.
    assert!(torus.fuse(OperationId(3), &rod).is_ok());
}
