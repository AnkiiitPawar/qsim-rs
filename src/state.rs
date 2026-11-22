use crate::complex::C64;
use crate::gates::{self, Gate};
use crate::rng::SplitMix64;

/// The full state of `n` qubits: 2^n complex amplitudes.
///
/// Qubit `q` is bit `q` of the basis-state index (little-endian), so the
/// basis state |q2 q1 q0⟩ has index q0 + 2 q1 + 4 q2.
#[derive(Clone, Debug, PartialEq)]
pub struct StateVector {
    n_qubits: usize,
    amps: Vec<C64>,
}

impl StateVector {
    /// |0…0⟩
    pub fn new(n_qubits: usize) -> Self {
        Self::basis(n_qubits, 0)
    }

    /// The computational basis state with the given index.
    pub fn basis(n_qubits: usize, index: usize) -> Self {
        assert!(n_qubits < usize::BITS as usize, "too many qubits");
        let mut amps = vec![C64::ZERO; 1 << n_qubits];
        amps[index] = C64::ONE;
        Self { n_qubits, amps }
    }

    pub fn from_amplitudes(amps: Vec<C64>) -> Self {
        assert!(
            amps.len().is_power_of_two(),
            "length must be a power of two"
        );
        Self {
            n_qubits: amps.len().trailing_zeros() as usize,
            amps,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    pub fn amplitudes(&self) -> &[C64] {
        &self.amps
    }

    /// Apply a 2×2 unitary to `target` when every qubit in `controls` is |1⟩.
    ///
    /// Each amplitude whose target bit is 0 is paired with the one whose target bit is 1,
    /// and the pair is multiplied by the matrix: O(2^n) work and no 2^n × 2^n matrix.
    pub fn apply_controlled(&mut self, controls: &[usize], target: usize, g: &Gate) -> &mut Self {
        assert!(target < self.n_qubits, "target out of range");
        let mut control_mask = 0usize;
        for &c in controls {
            assert!(c < self.n_qubits && c != target, "bad control qubit");
            control_mask |= 1 << c;
        }
        let stride = 1usize << target;
        let dim = self.amps.len();
        let mut base = 0;
        while base < dim {
            for i in base..base + stride {
                if i & control_mask != control_mask {
                    continue;
                }
                let j = i | stride;
                let (a, b) = (self.amps[i], self.amps[j]);
                self.amps[i] = g[0][0] * a + g[0][1] * b;
                self.amps[j] = g[1][0] * a + g[1][1] * b;
            }
            base += 2 * stride;
        }
        self
    }

    pub fn apply(&mut self, target: usize, g: &Gate) -> &mut Self {
        self.apply_controlled(&[], target, g)
    }

    pub fn h(&mut self, q: usize) -> &mut Self {
        self.apply(q, &gates::H)
    }
    pub fn x(&mut self, q: usize) -> &mut Self {
        self.apply(q, &gates::X)
    }
    pub fn y(&mut self, q: usize) -> &mut Self {
        self.apply(q, &gates::Y)
    }
    pub fn z(&mut self, q: usize) -> &mut Self {
        self.apply(q, &gates::Z)
    }
    pub fn cnot(&mut self, control: usize, target: usize) -> &mut Self {
        self.apply_controlled(&[control], target, &gates::X)
    }
    pub fn cz(&mut self, control: usize, target: usize) -> &mut Self {
        self.apply_controlled(&[control], target, &gates::Z)
    }
    /// Controlled phase: multiplies |11⟩ on (control, target) by e^{iθ}.
    pub fn cphase(&mut self, control: usize, target: usize, theta: f64) -> &mut Self {
        self.apply_controlled(&[control], target, &gates::phase(theta))
    }
    pub fn toffoli(&mut self, c1: usize, c2: usize, target: usize) -> &mut Self {
        self.apply_controlled(&[c1, c2], target, &gates::X)
    }

    pub fn swap(&mut self, a: usize, b: usize) -> &mut Self {
        assert!(a < self.n_qubits && b < self.n_qubits);
        if a != b {
            let (ma, mb) = (1usize << a, 1usize << b);
            for i in 0..self.amps.len() {
                if i & ma != 0 && i & mb == 0 {
                    self.amps.swap(i, i ^ ma ^ mb);
                }
            }
        }
        self
    }

    /// Multiply one basis state's amplitude by -1 (a phase oracle).
    pub fn flip_phase(&mut self, index: usize) -> &mut Self {
        self.amps[index] = -self.amps[index];
        self
    }

    /// 2|0⟩⟨0| - I: keep |0…0⟩, negate everything else.
    pub fn reflect_about_zero(&mut self) -> &mut Self {
        for a in self.amps.iter_mut().skip(1) {
            *a = -*a;
        }
        self
    }

    pub fn probabilities(&self) -> Vec<f64> {
        self.amps.iter().map(|a| a.norm_sqr()).collect()
    }

    pub fn norm_sqr(&self) -> f64 {
        self.amps.iter().map(|a| a.norm_sqr()).sum()
    }

    /// ⟨self|other⟩
    pub fn inner(&self, other: &Self) -> C64 {
        let mut acc = C64::ZERO;
        for (a, b) in self.amps.iter().zip(&other.amps) {
            acc += a.conj() * *b;
        }
        acc
    }

    /// Sample a basis state from the Born-rule distribution, without collapsing.
    pub fn sample(&self, rng: &mut SplitMix64) -> usize {
        let mut u = rng.next_f64();
        for (i, a) in self.amps.iter().enumerate() {
            u -= a.norm_sqr();
            if u < 0.0 {
                return i;
            }
        }
        self.amps.len() - 1
    }

    /// Measure one qubit: returns the outcome and collapses (renormalises) the state.
    pub fn measure(&mut self, q: usize, rng: &mut SplitMix64) -> bool {
        let mask = 1usize << q;
        let p1: f64 = self
            .amps
            .iter()
            .enumerate()
            .filter(|(i, _)| i & mask != 0)
            .map(|(_, a)| a.norm_sqr())
            .sum();
        let outcome = rng.next_f64() < p1;
        let keep_norm = if outcome { p1 } else { 1.0 - p1 }.sqrt();
        for (i, a) in self.amps.iter_mut().enumerate() {
            *a = if (i & mask != 0) == outcome {
                *a * (1.0 / keep_norm)
            } else {
                C64::ZERO
            };
        }
        outcome
    }
}
