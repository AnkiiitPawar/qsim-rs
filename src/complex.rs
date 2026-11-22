use std::ops::{Add, AddAssign, Mul, Neg, Sub};

/// A complex number with f64 parts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub const ZERO: C64 = C64::new(0.0, 0.0);
    pub const ONE: C64 = C64::new(1.0, 0.0);
    pub const I: C64 = C64::new(0.0, 1.0);

    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    /// r e^{iθ}
    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.norm_sqr().sqrt()
    }
}

impl Add for C64 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}

impl AddAssign for C64 {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Sub for C64 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for C64 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl Mul<f64> for C64 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}

impl Neg for C64 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
