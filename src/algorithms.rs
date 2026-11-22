//! Standard circuits built from the gate primitives.

use crate::state::StateVector;
use std::f64::consts::PI;

/// (|00⟩ + |11⟩) / √2
pub fn bell() -> StateVector {
    let mut s = StateVector::new(2);
    s.h(0).cnot(0, 1);
    s
}

/// (|0…0⟩ + |1…1⟩) / √2 on n qubits.
pub fn ghz(n: usize) -> StateVector {
    let mut s = StateVector::new(n);
    s.h(0);
    for q in 1..n {
        s.cnot(q - 1, q);
    }
    s
}

/// Quantum Fourier transform on all qubits:
/// |x⟩ → (1/√N) Σ_y e^{2πi xy/N} |y⟩.
///
/// Uses n Hadamards, n(n-1)/2 controlled phases and ⌊n/2⌋ swaps: O(n²) gates,
/// compared with the O(N log N) classical FFT on N = 2^n amplitudes.
pub fn qft(s: &mut StateVector) {
    let n = s.n_qubits();
    for target in (0..n).rev() {
        s.h(target);
        for control in (0..target).rev() {
            s.cphase(control, target, PI / (1u64 << (target - control)) as f64);
        }
    }
    for q in 0..n / 2 {
        s.swap(q, n - 1 - q);
    }
}

/// Exact inverse of [`qft`]: the same gates in reverse order with conjugated phases.
pub fn inverse_qft(s: &mut StateVector) {
    let n = s.n_qubits();
    for q in 0..n / 2 {
        s.swap(q, n - 1 - q);
    }
    for target in 0..n {
        for control in 0..target {
            s.cphase(control, target, -PI / (1u64 << (target - control)) as f64);
        }
        s.h(target);
    }
}

/// Number of Grover iterations that maximises the success probability:
/// ⌊π / (4θ)⌋ with sin θ = 1/√N.
pub fn grover_iterations(n_qubits: usize) -> usize {
    let theta = (1.0 / ((1u64 << n_qubits) as f64).sqrt()).asin();
    (PI / (4.0 * theta)).floor() as usize
}

/// One Grover iteration: phase oracle for `marked`, then the diffusion
/// operator H^⊗n (2|0⟩⟨0| - I) H^⊗n = 2|s⟩⟨s| - I.
pub fn grover_iteration(s: &mut StateVector, marked: usize) {
    let n = s.n_qubits();
    s.flip_phase(marked);
    for q in 0..n {
        s.h(q);
    }
    s.reflect_about_zero();
    for q in 0..n {
        s.h(q);
    }
}

/// Grover search for one marked item among N = 2^n, run for `iterations` steps.
pub fn grover(n_qubits: usize, marked: usize, iterations: usize) -> StateVector {
    let mut s = StateVector::new(n_qubits);
    for q in 0..n_qubits {
        s.h(q);
    }
    for _ in 0..iterations {
        grover_iteration(&mut s, marked);
    }
    s
}
