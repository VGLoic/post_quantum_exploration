/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use super::polynomial::Polynomial;
use crate::{modular_arithmetic, primefield::PrimeFieldElement};

const N: u32 = 1_000_000_007;

type Polynomial1B7 = Polynomial<N>;

fn generate_proof() {
    let constraint_polynomial =
        Polynomial1B7::interpolate_from_roots((0..=9).map(PrimeFieldElement::from).collect());
    let z_polynomial = Polynomial1B7::interpolate_from_roots(
        (0..=1_000_000_000).map(PrimeFieldElement::from).collect(),
    );
}

#[cfg(test)]
mod test {
    use super::*;
}
