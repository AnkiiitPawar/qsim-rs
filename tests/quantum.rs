use qsim::algorithms::{bell, ghz, grover, grover_iterations, inverse_qft, qft};
use qsim::{gates, SplitMix64, StateVector, C64};
use std::f64::consts::PI;

const TOL: f64 = 1e-12;

fn random_state(n: usize, seed: u64) -> StateVector {
    let mut rng = SplitMix64::new(seed);
    let amps: Vec<C64> = (0..1 << n)
        .map(|_| C64::new(rng.next_f64() - 0.5, rng.next_f64() - 0.5))
        .collect();
    let norm = amps.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
    StateVector::from_amplitudes(amps.into_iter().map(|a| a * (1.0 / norm)).collect())
}

fn assert_close(a: &StateVector, b: &StateVector) {
    for (x, y) in a.amplitudes().iter().zip(b.amplitudes()) {
        assert!((*x - *y).abs() < 1e-10, "{x:?} vs {y:?}");
    }
}

#[test]
fn hadamard_is_its_own_inverse() {
    let original = random_state(5, 1);
    let mut s = original.clone();
    s.h(3).h(3);
    assert_close(&s, &original);
}

#[test]
fn pauli_algebra_xyz_is_i_times_identity() {
    // XYZ = iI, so applying Z, then Y, then X multiplies every amplitude by i.
    let original = random_state(3, 2);
    let mut s = original.clone();
    s.z(1).y(1).x(1);
    for (x, y) in s.amplitudes().iter().zip(original.amplitudes()) {
        assert!((*x - C64::I * *y).abs() < TOL);
    }
}

#[test]
fn bell_and_ghz_states_have_the_right_correlations() {
    let p = bell().probabilities();
    assert!((p[0] - 0.5).abs() < TOL && (p[3] - 0.5).abs() < TOL);
    assert!(p[1] < TOL && p[2] < TOL);

    let p = ghz(6).probabilities();
    assert!((p[0] - 0.5).abs() < TOL && (p[63] - 0.5).abs() < TOL);
    assert!((p.iter().sum::<f64>() - 1.0).abs() < TOL);
}

#[test]
fn toffoli_truth_table() {
    for input in 0..8usize {
        let mut s = StateVector::basis(3, input);
        s.toffoli(0, 1, 2);
        let expected = if input & 0b011 == 0b011 {
            input ^ 0b100
        } else {
            input
        };
        assert!(
            (s.probabilities()[expected] - 1.0).abs() < TOL,
            "input {input:03b}"
        );
    }
}

#[test]
fn qft_matches_the_discrete_fourier_transform() {
    let n = 4;
    let dim = 1usize << n;
    for x in 0..dim {
        let mut s = StateVector::basis(n, x);
        qft(&mut s);
        for (y, amp) in s.amplitudes().iter().enumerate() {
            let expected = C64::from_polar(
                1.0 / (dim as f64).sqrt(),
                2.0 * PI * (x * y) as f64 / dim as f64,
            );
            assert!(
                (*amp - expected).abs() < 1e-12,
                "x={x} y={y}: {amp:?} vs {expected:?}"
            );
        }
    }
}

#[test]
fn inverse_qft_undoes_qft() {
    let original = random_state(7, 3);
    let mut s = original.clone();
    qft(&mut s);
    inverse_qft(&mut s);
    assert_close(&s, &original);
}

#[test]
fn grover_matches_the_analytic_success_probability() {
    // After r iterations the marked item has probability sin²((2r + 1)θ), sin θ = 1/√N.
    let n = 6;
    let marked = 0b101101;
    let theta = (1.0 / ((1 << n) as f64).sqrt()).asin();
    for r in 0..12 {
        let p = grover(n, marked, r).probabilities()[marked];
        let expected = ((2 * r + 1) as f64 * theta).sin().powi(2);
        assert!((p - expected).abs() < 1e-12, "r={r}: {p} vs {expected}");
    }
    let best = grover_iterations(n);
    assert_eq!(best, 6);
    assert!(grover(n, marked, best).probabilities()[marked] > 0.99);
}

#[test]
fn random_circuits_preserve_the_norm() {
    let mut rng = SplitMix64::new(4);
    let n = 8;
    let mut s = StateVector::new(n);
    for _ in 0..500 {
        let q = (rng.next_u64() % n as u64) as usize;
        let c = (q + 1 + (rng.next_u64() % (n as u64 - 1)) as usize) % n;
        let angle = rng.next_f64() * 2.0 * PI;
        match rng.next_u64() % 6 {
            0 => s.h(q),
            1 => s.apply(q, &gates::rx(angle)),
            2 => s.apply(q, &gates::ry(angle)),
            3 => s.apply(q, &gates::rz(angle)),
            4 => s.cnot(c, q),
            _ => s.cphase(c, q, angle),
        };
    }
    assert!((s.norm_sqr() - 1.0).abs() < 1e-12);
}

#[test]
fn measurement_statistics_follow_the_born_rule() {
    let mut rng = SplitMix64::new(5);
    let s = bell();
    let shots = 20_000;
    let mut counts = [0usize; 4];
    for _ in 0..shots {
        counts[s.sample(&mut rng)] += 1;
    }
    assert_eq!(counts[1] + counts[2], 0);
    // binomial standard deviation sqrt(N p (1-p)) ≈ 71
    assert!((counts[0] as f64 - shots as f64 / 2.0).abs() < 4.0 * 71.0);
}

#[test]
fn measuring_one_half_of_a_bell_pair_fixes_the_other() {
    let mut rng = SplitMix64::new(6);
    for _ in 0..100 {
        let mut s = bell();
        let first = s.measure(0, &mut rng);
        let second = s.measure(1, &mut rng);
        assert_eq!(first, second);
        assert!((s.norm_sqr() - 1.0).abs() < TOL);
    }
}
