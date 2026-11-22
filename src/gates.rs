//! Single-qubit gates as 2×2 unitary matrices, written row-major as [[m00, m01], [m10, m11]].

use crate::complex::C64;
use std::f64::consts::FRAC_1_SQRT_2;

pub type Gate = [[C64; 2]; 2];

const fn re(x: f64) -> C64 {
    C64::new(x, 0.0)
}

pub const X: Gate = [[re(0.0), re(1.0)], [re(1.0), re(0.0)]];
pub const Y: Gate = [
    [re(0.0), C64::new(0.0, -1.0)],
    [C64::new(0.0, 1.0), re(0.0)],
];
pub const Z: Gate = [[re(1.0), re(0.0)], [re(0.0), re(-1.0)]];
pub const H: Gate = [
    [re(FRAC_1_SQRT_2), re(FRAC_1_SQRT_2)],
    [re(FRAC_1_SQRT_2), re(-FRAC_1_SQRT_2)],
];
pub const S: Gate = [[re(1.0), re(0.0)], [re(0.0), C64::new(0.0, 1.0)]];

/// T = diag(1, e^{iπ/4})
pub fn t() -> Gate {
    phase(std::f64::consts::FRAC_PI_4)
}

/// diag(1, e^{iθ})
pub fn phase(theta: f64) -> Gate {
    [[re(1.0), re(0.0)], [re(0.0), C64::from_polar(1.0, theta)]]
}

/// exp(-iθX/2)
pub fn rx(theta: f64) -> Gate {
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    [[re(c), C64::new(0.0, -s)], [C64::new(0.0, -s), re(c)]]
}

/// exp(-iθY/2)
pub fn ry(theta: f64) -> Gate {
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    [[re(c), re(-s)], [re(s), re(c)]]
}

/// exp(-iθZ/2)
pub fn rz(theta: f64) -> Gate {
    [
        [C64::from_polar(1.0, -theta / 2.0), re(0.0)],
        [re(0.0), C64::from_polar(1.0, theta / 2.0)],
    ]
}
