/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use crate::primefield::PrimeFieldElement;

#[derive(Clone, Debug)]
pub struct StarkLeaf<const N: u32> {
    pub cp: PrimeFieldElement<N>,
    pub d: PrimeFieldElement<N>,
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
    use crate::{
        merkletree_v2::{MerkleTreeV2, ValueWithProof, verify_proof},
        stark::polynomial::Polynomial,
    };

    const N: u32 = 1_000_000_007;
    const P_DEGREE: u32 = 1_000; // 0 <= P(x) <= 9 for 0 <= x <= 1_000 => 1_000 constraints, we can easily satisfy it with a 1_000 degree polynomial
    const TOTAL_POINTS: u32 = 1_000_000;

    type FieldElement1B7 = PrimeFieldElement<N>;
    type Polynomial1B7 = Polynomial<N>;

    fn generate_commitments_tree() -> Result<MerkleTreeV2<StarkLeaf<N>>, anyhow::Error> {
        let mut points = Vec::with_capacity(1 + P_DEGREE as usize);
        let mut values = Vec::with_capacity(1 + P_DEGREE as usize);
        let mut z_roots = Vec::with_capacity(1 + P_DEGREE as usize);

        for i in 0_u32..=P_DEGREE {
            let i_as_field_element = FieldElement1B7::from(i);
            let y = rand::random_range(0u32..10u32);
            values.push(FieldElement1B7::from(y));
            points.push(i_as_field_element);
            z_roots.push(i_as_field_element);
        }

        let p = Polynomial1B7::interpolate_from_coordinates(points, values).ok_or(anyhow!(
            "Unable to interpolate P polynomial from coordinates"
        ))?;

        if p.degree() != P_DEGREE as usize {
            return Err(anyhow!("Invalid degree, got {}", p.degree()));
        }
        let constraint_polynomial =
            Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
        let z_polynomial = Polynomial1B7::interpolate_from_roots(z_roots);

        let mut values: Vec<StarkLeaf<N>> = Vec::with_capacity(TOTAL_POINTS as usize);
        for i in 0u32..TOTAL_POINTS {
            let i_as_field_element = FieldElement1B7::from(i);

            let p_eval = p.evaluate(i_as_field_element);
            let cp_eval = constraint_polynomial.evaluate(p_eval);
            let d_eval = if i <= P_DEGREE {
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

        MerkleTreeV2::new(20, &values)
    }

    fn point_to_selector(point: PrimeFieldElement<N>) -> Vec<bool> {
        let mut a = point.inner();
        let mut selector = vec![];

        while a != 0 {
            if a % 2 == 0 {
                selector.push(false);
            } else {
                selector.push(true);
            }
            a /= 2;
        }

        if selector.len() < 20 {
            selector.resize(20, false);
        }

        selector.reverse();

        selector
    }

    #[test]
    fn test_selector() {
        assert_eq!(point_to_selector(0.into()), vec![false; 20]);
        assert_eq!(point_to_selector(1_048_575.into()), vec![true; 20]);
        let mut expected = vec![false; 20];
        expected[19] = true;
        assert_eq!(point_to_selector(1.into()), expected);
        let mut expected = vec![false; 20];
        expected[18] = true;
        assert_eq!(point_to_selector(2.into()), expected);
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        // ########################
        // ###### Alice part ######
        // ########################

        let commitments_tree = generate_commitments_tree().unwrap();

        // Bob asks 16 points
        let mut proofs: Vec<(PrimeFieldElement<N>, ValueWithProof<StarkLeaf<N>>)> = vec![];
        for _ in 0..16 {
            let point_u32: u32 = rand::random_range(0..TOTAL_POINTS);
            let point = PrimeFieldElement::<N>::from(point_u32);
            proofs.push((
                point,
                commitments_tree.get(&point_to_selector(point)).unwrap(),
            ));
        }

        // ######################
        // ###### Bob part ######
        // ######################

        for (point, value_with_proof) in proofs {
            if point.inner() <= P_DEGREE {
                assert_eq!(value_with_proof.value.cp, 0.into());
            }
            let formatted_proof: Vec<_> = value_with_proof
                .proof
                .into_iter()
                .map(|p| (p.0, p.1.to_owned()))
                .collect();
            verify_proof(
                commitments_tree.root_hash(),
                value_with_proof.value,
                &formatted_proof,
            )
            .expect("invalid proof");
        }
    }
}
