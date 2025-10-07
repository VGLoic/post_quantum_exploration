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
    use sha3::{Digest, Sha3_256};

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

    fn generate_commitments_tree(
        p: &Polynomial1B7,
        constraint_polynomial: &Polynomial1B7,
        z_polynomial: &Polynomial1B7,
    ) -> Result<MerkleTreeV2<StarkLeaf<N>>, anyhow::Error> {
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

    fn derive_proofs<'a>(
        tree: &'a MerkleTreeV2<StarkLeaf<N>>,
    ) -> Result<Vec<(PrimeFieldElement<N>, ValueWithProof<'a, StarkLeaf<N>>)>, anyhow::Error> {
        let points = root_to_points(tree.root_hash());
        let mut proofs: Vec<(PrimeFieldElement<N>, ValueWithProof<StarkLeaf<N>>)> = vec![];
        for point in points {
            let selector = tree.get(&point_to_selector(point))?;
            proofs.push((point, selector))
        }
        Ok(proofs)
    }

    fn root_to_points(root: &[u8; 32]) -> Vec<PrimeFieldElement<N>> {
        let mut points = Vec::with_capacity(16);
        for chunk in root.chunks(2) {
            let digest = Sha3_256::digest(chunk);
            let mut le_bytes = [0u8; 4];
            le_bytes.clone_from_slice(&digest[0..4]);
            let point_as_u32 = u32::from_le_bytes(le_bytes) % TOTAL_POINTS;
            points.push(PrimeFieldElement::<N>::from(point_as_u32));
        }
        points
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
        // ###################
        // ###### Setup ######
        // ###################

        let points: Vec<FieldElement1B7> = (0..(1 + P_DEGREE)).map(FieldElement1B7::from).collect();
        let constraint_polynomial =
            Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
        let z_polynomial = Polynomial1B7::interpolate_from_roots(points.clone());

        // ########################
        // ###### Alice part ######
        // ########################

        let p_values = (0..(1 + P_DEGREE))
            .map(|_| FieldElement1B7::from(rand::random_range(0u32..10u32)))
            .collect();
        let p = Polynomial1B7::interpolate_from_coordinates(points, p_values)
            .ok_or(anyhow!(
                "Unable to interpolate P polynomial from coordinates"
            ))
            .unwrap();

        if p.degree() != P_DEGREE as usize {
            panic!("Invalid degree, got {}", p.degree());
        }

        let commitments_tree =
            generate_commitments_tree(&p, &constraint_polynomial, &z_polynomial).unwrap();
        let proofs = derive_proofs(&commitments_tree).unwrap();

        // ######################
        // ###### Bob part ######
        // ######################

        let expected_points = root_to_points(commitments_tree.root_hash());

        for ((point, value_with_proof), expected_point) in proofs.iter().zip(&expected_points) {
            assert_eq!(point, expected_point);
            if point.inner() <= P_DEGREE {
                assert_eq!(value_with_proof.value.cp, 0.into());
            } else {
                let z_eval = z_polynomial.evaluate(*point);
                assert_eq!(
                    z_eval.mul(&value_with_proof.value.d),
                    value_with_proof.value.cp
                );
            }
            let formatted_proof: Vec<_> = value_with_proof
                .proof
                .iter()
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
