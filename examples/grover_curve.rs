//! Success probability of Grover search vs number of iterations (CSV on stdout),
//! plus the final probability distribution after the optimal number of steps.
//! Run with: cargo run --release --example grover_curve > docs/grover.csv

use qsim::algorithms::{grover_iteration, grover_iterations};
use qsim::StateVector;

fn main() {
    let n = 8;
    let marked = 0b1011_0110;
    let best = grover_iterations(n);
    let theta = (1.0 / ((1 << n) as f64).sqrt()).asin();

    let mut s = StateVector::new(n);
    for q in 0..n {
        s.h(q);
    }
    println!("kind,x,simulated,theory");
    for r in 0..=3 * best {
        let p = s.probabilities()[marked];
        let theory = ((2 * r + 1) as f64 * theta).sin().powi(2);
        println!("curve,{r},{p:.12},{theory:.12}");
        if r == best {
            for (k, pk) in s.probabilities().iter().enumerate() {
                println!("dist,{k},{pk:.12},");
            }
        }
        grover_iteration(&mut s, marked);
    }
    eprintln!(
        "N = {}, optimal iterations = {best}, marked = {marked}",
        1 << n
    );
}
