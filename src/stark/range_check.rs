/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use crate::primefield::PrimeFieldElement;

#[derive(Clone, Debug)]
pub struct StarkLeaf<const N: u32> {
    cp: PrimeFieldElement<N>,
    d: PrimeFieldElement<N>,
    rep: [u8; 8],
}

impl<const N: u32> StarkLeaf<N> {
    pub fn new(cp: PrimeFieldElement<N>, d: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 8];
        rep[0..4].clone_from_slice(&cp.inner().to_le_bytes());
        rep[4..].clone_from_slice(&d.inner().to_le_bytes());
        Self { cp, d, rep }
    }
}

impl<const N: u32> AsRef<[u8]> for StarkLeaf<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u32> Default for StarkLeaf<N> {
    fn default() -> Self {
        Self {
            cp: 0.into(),
            d: 0.into(),
            rep: [0u8; 8],
        }
    }
}

#[cfg(test)]
mod test {
    use anyhow::anyhow;

    use super::*;
    use crate::{merkletree_v2::MerkleTreeV2, stark::polynomial::Polynomial};

    const N: u32 = 1_000_000_007;
    const P_DEGREE: usize = 100;
    const Z_DEGREE: usize = 1_000;
    const TOTAL_POINTS: u32 = 1_000_000;

    type FieldElement1B7 = PrimeFieldElement<N>;
    type Polynomial1B7 = Polynomial<N>;

    fn generate_proof() -> Result<u32, anyhow::Error> {
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

        let p = Polynomial1B7::interpolate_from_coordinates(points, values).ok_or(anyhow!(
            "Unable to interpolate P polynomial from coordinates"
        ))?;
        let constraint_polynomial =
            Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
        let z_polynomial = Polynomial1B7::interpolate_from_roots(z_roots);

        let mut values: Vec<StarkLeaf<N>> = Vec::with_capacity(1 + Z_DEGREE);
        for i in 0u32..TOTAL_POINTS {
            let i_as_field_element = FieldElement1B7::from(i);

            let p_eval = p.evaluate(i_as_field_element);
            let cp_eval = constraint_polynomial.evaluate(p_eval);
            let d_eval = if (i as usize) < Z_DEGREE {
                0.into()
            } else {
                let z_eval = z_polynomial.evaluate(i_as_field_element);
                cp_eval.mul(
                    &z_eval
                        .inv()
                        .ok_or(anyhow!("unable to find inverse of Z evaluation"))?,
                )
            };

            values.push(StarkLeaf::new(cp_eval, d_eval));
        }

        let _tree = MerkleTreeV2::new(20, &values)?;

        Ok(3)
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        let _proof = generate_proof().unwrap();
    }
}
