/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use super::polynomial::Polynomial;
use crate::{modular_arithmetic, primefield::PrimeFieldElement};

const N: u32 = 1_000_000_007;

type Polynomial1B7 = Polynomial<N>;
type FieldElement1B7 = PrimeFieldElement<N>;

#[derive(Clone, Debug)]
struct StarkLeaf {
    cp: FieldElement1B7,
    d: FieldElement1B7,
    rep: [u8; 64],
}

impl AsRef<[u8]> for StarkLeaf {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl Default for StarkLeaf {
    fn default() -> Self {
        Self {
            cp: 0.into(),
            d: 0.into(),
            rep: [0u8; 64],
        }
    }
}

const P_DEGREE: usize = 10_000;
const Z_DEGREE: usize = 1_000_000_000;

pub fn generate_proof() -> Option<u32> {
    let mut points = Vec::with_capacity(P_DEGREE);
    let mut values = Vec::with_capacity(P_DEGREE);
    let mut z_roots = Vec::with_capacity(Z_DEGREE);

    for i in 0_u32..(Z_DEGREE as u32) {
        let i_as_field_element = FieldElement1B7::from(i);
        if (i as usize) < P_DEGREE {
            points.push(i_as_field_element);
            let y = rand::random_range(5u32..2_000u32);
            values.push(FieldElement1B7::from(y));
        }
        z_roots.push(i_as_field_element);
    }

    let p = Polynomial1B7::interpolate_from_coordinates(points, values)?;
    // let constraint_polynomial =
    //     Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
    // let z_polynomial = Polynomial1B7::interpolate_from_roots(z_roots);

    Some(3)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        let proof = generate_proof();
        assert!(proof.is_some());
    }
}
