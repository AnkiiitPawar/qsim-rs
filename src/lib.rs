//! A small state-vector quantum circuit simulator.
//!
//! ```
//! use qsim::{algorithms, StateVector};
//!
//! let bell = algorithms::bell();
//! let p = bell.probabilities();
//! assert!((p[0b00] - 0.5).abs() < 1e-12 && (p[0b11] - 0.5).abs() < 1e-12);
//!
//! let mut s = StateVector::new(3);
//! s.h(0).cnot(0, 1).toffoli(0, 1, 2); // chain gates
//! ```

pub mod algorithms;
pub mod complex;
pub mod gates;
pub mod rng;
pub mod state;

pub use complex::C64;
pub use rng::SplitMix64;
pub use state::StateVector;
