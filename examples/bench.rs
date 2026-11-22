//! Time per gate as the register grows: every extra qubit doubles the state vector.
//! Run with: cargo run --release --example bench

use qsim::StateVector;
use std::time::Instant;

fn main() {
    println!(
        "{:>7} {:>14} {:>14} {:>16}",
        "qubits", "amplitudes", "memory", "ns per H gate"
    );
    for n in (12..=24).step_by(2) {
        let mut s = StateVector::new(n);
        let reps = (1 << 22) / (1 << n.min(22)) + 1;
        let start = Instant::now();
        for _ in 0..reps {
            for q in 0..n {
                s.h(q);
            }
        }
        let per_gate = start.elapsed().as_nanos() as f64 / (reps * n) as f64;
        let mib = (16u64 << n) as f64 / (1 << 20) as f64;
        println!(
            "{n:>7} {:>14} {:>11.1} MiB {per_gate:>16.0}",
            1u64 << n,
            mib
        );
        assert!((s.norm_sqr() - 1.0).abs() < 1e-9);
    }
}
