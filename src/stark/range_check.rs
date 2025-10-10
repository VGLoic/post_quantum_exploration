/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use crate::primefield::PrimeFieldElement;

#[derive(Clone, Debug)]
pub struct OriginalEvaluation<const N: u32> {
    pub p: PrimeFieldElement<N>,
    pub d: PrimeFieldElement<N>,
    rep: [u8; 8],
}

impl<const N: u32> OriginalEvaluation<N> {
    pub fn new(p: PrimeFieldElement<N>, d: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 8];
        rep[0..4].clone_from_slice(&p.inner().to_le_bytes());
        rep[4..].clone_from_slice(&d.inner().to_le_bytes());
        Self { p, d, rep }
    }
}

impl<const N: u32> AsRef<[u8]> for OriginalEvaluation<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u32> Default for OriginalEvaluation<N> {
    fn default() -> Self {
        Self {
            p: 0.into(),
            d: 0.into(),
            rep: [0u8; 8],
        }
    }
}

#[derive(Clone, Debug)]
pub struct ColumnEvaluation<const N: u32> {
    pub v: PrimeFieldElement<N>,
    rep: [u8; 4],
}

impl<const N: u32> ColumnEvaluation<N> {
    pub fn new(v: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 4];
        rep[0..4].clone_from_slice(&v.inner().to_le_bytes());
        Self { v, rep }
    }
}

impl<const N: u32> AsRef<[u8]> for ColumnEvaluation<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u32> Default for ColumnEvaluation<N> {
    fn default() -> Self {
        Self {
            v: 0.into(),
            rep: [0u8; 4],
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use anyhow::anyhow;
    use sha3::{Digest, Sha3_256};

    use super::*;
    use crate::{
        merkletree_v2::{MerkleProof, MerkleTreeV2, verify_proof},
        stark::polynomial::Polynomial,
    };

    const N: u32 = 1_048_589; // Chosen because: (p - 1) / 4 = 262147 and p > 2**20 (depth of the tree)
    const MINUS_ONE_SQRT_0: u32 = 38993; // sqrt(-1) (mod N)
    const MINUS_ONE_SQRT_1: u32 = 1009596; // sqrt(-1) (mod N)
    const P_MAX_DEGREE: u32 = 1_000; // 0 <= P(x) <= 9 for 0 <= x <= 999 => 1_000 constraints. One can always use Lagrange to build a valid 1_000 degree polynomial
    const TOTAL_POINTS: u32 = 1_048_576; // We take this in order to completely fill a Merkle tree of depth 20 while being lower than N

    type FieldElementRG = PrimeFieldElement<N>;
    type PolynomialRG = Polynomial<N>;

    fn point_to_selector(point: PrimeFieldElement<N>, depth: usize) -> Vec<bool> {
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

        if selector.len() < depth {
            selector.resize(depth, false);
        }

        selector.reverse();

        selector
    }

    fn pseudo_random_select_column_point(root: &[u8; 32], modulus: u32) -> PrimeFieldElement<N> {
        let chunk = &root[0..2];
        let digest = Sha3_256::digest(chunk);
        let mut le_bytes = [0u8; 4];
        le_bytes.clone_from_slice(&digest[0..4]);
        let point_as_u32 = u32::from_le_bytes(le_bytes) % modulus;
        point_as_u32.into()
    }

    fn pseudo_random_select_column_points(
        root: &[u8; 32],
        root_4_groups: &[FourRoot],
    ) -> Vec<(u32, FieldElementRG)> {
        let mut points = vec![];
        let mut selected_indices = HashSet::new();
        let mut counter = 0;
        let mut seed: [u8; 32] = *root;
        while points.len() < 300 {
            seed = Sha3_256::digest(seed).into(); // 32 bytes
            for chunk in seed.chunks(4) {
                let mut le_bytes = [0u8; 4];
                le_bytes.clone_from_slice(chunk);

                let root_4_group_index =
                    u32::from_le_bytes(le_bytes) % (root_4_groups.len() as u32);
                if !selected_indices.contains(&root_4_group_index) {
                    let root_4_group = &root_4_groups[root_4_group_index as usize];
                    if root_4_group.iter().all(|el| el.inner() < TOTAL_POINTS) {
                        points.push((root_4_group_index, root_4_group.first()));
                        selected_indices.insert(root_4_group_index);
                    }
                }
            }
            counter += 1;

            if counter > 10_000 {
                panic!("oula");
            }
        }
        points.resize(300, (0, 0.into()));
        points
    }

    struct FourRoot([FieldElementRG; 4]);
    impl FourRoot {
        fn first(&self) -> FieldElementRG {
            self.0[0]
        }
        fn iter(&self) -> impl Iterator<Item = &FieldElementRG> {
            self.0.iter()
        }
        fn four_exp(&self) -> FieldElementRG {
            self.0[0].exp(4)
        }
    }

    fn group_roots_4(max: u32) -> Vec<FourRoot> {
        let zero = PrimeFieldElement::from(0);
        let mut elements = vec![FourRoot([zero; 4])];

        let mut covered_elements = HashSet::new();

        let minus_one_sqrt_0 = FieldElementRG::from(MINUS_ONE_SQRT_0);
        let minus_one_sqrt_1 = FieldElementRG::from(MINUS_ONE_SQRT_1);

        for i in 1..max {
            let i_as_field_element = FieldElementRG::from(i);
            if !covered_elements.contains(&i_as_field_element) {
                let root_1 = i_as_field_element.mul(&minus_one_sqrt_0);
                let root_2 = i_as_field_element.mul(&minus_one_sqrt_1);
                let root_3 = i_as_field_element.neg();
                covered_elements.insert(root_1);
                covered_elements.insert(root_2);
                covered_elements.insert(root_3);
                elements.push(FourRoot([i_as_field_element, root_1, root_2, root_3]));
            }
        }

        elements
    }

    #[test]
    fn test_selector() {
        assert_eq!(point_to_selector(0.into(), 20), vec![false; 20]);
        assert_eq!(point_to_selector(1_048_575.into(), 20), vec![true; 20]);
        let mut expected = vec![false; 20];
        expected[19] = true;
        assert_eq!(point_to_selector(1.into(), 20), expected);
        let mut expected = vec![false; 20];
        expected[18] = true;
        assert_eq!(point_to_selector(2.into(), 20), expected);
    }

    #[test]
    fn test_small_image() {
        let root_groups = group_roots_4(N);
        assert_eq!(root_groups.len() as u32, (N - 1) / 4 + 1);
        for root_group in root_groups {
            let expected = root_group.four_exp();
            for v in root_group.iter().skip(1).map(|el| el.exp(4)) {
                assert_eq!(v, expected);
            }
        }
    }

    #[derive(Debug)]
    struct OriginalCommitment<const N: u32> {
        point: PrimeFieldElement<N>,
        evaluation: OriginalEvaluation<N>,
        proof: MerkleProof,
    }

    #[derive(Debug)]
    struct ColumnCommitment<const N: u32> {
        point: PrimeFieldElement<N>,
        evaluation: ColumnEvaluation<N>,
        proof: MerkleProof,
    }

    type FourRootOriginalCommitments<const N: u32> = [OriginalCommitment<N>; 4];

    struct StarkProof<const N: u32> {
        pub original_commitment_root: [u8; 32],
        pub original: Vec<FourRootOriginalCommitments<N>>,
        pub column_root: [u8; 32],
        pub column: Vec<ColumnCommitment<N>>,
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        // ###################
        // ###### Setup ######
        // ###################

        let constraint_polynomial =
            PolynomialRG::interpolate_from_roots((0..=9).map(FieldElementRG::from).collect()); // x(x - 1)...(x - 9)

        // ########################
        // ###### Alice part ######
        // ########################

        // We define a polynomial P with degree less than P_MAX_DEGREE and verifying 0 <= P(x) <= 9 for x between 0 and <P_MAX_DEGREE
        // P will be 0 at 0 <= x < P_MAX_DEGREE - 1, it will gives a (P_MAX_DEGREE - 1) polynomial that we're gonna scale down with the evaluation at x = P_MAX_DEGREE - 1
        let mut p = PolynomialRG::interpolate_from_roots(
            (0..(P_MAX_DEGREE - 1)).map(FieldElementRG::from).collect(),
        );
        p = p.mul_by_scalar(p.evaluate((P_MAX_DEGREE - 1).into()).inv().unwrap());
        if p.degree() >= P_MAX_DEGREE as usize {
            panic!("Invalid degree, got {}", p.degree());
        }

        // We evaluate the polynomial over 0..TOTAL_POINTS
        let mut original_evaluations: Vec<OriginalEvaluation<N>> =
            Vec::with_capacity(TOTAL_POINTS as usize);
        for i in 0u32..TOTAL_POINTS {
            let i_as_field_element = FieldElementRG::from(i);

            let p_eval = p.evaluate(i_as_field_element);
            let cp_eval = constraint_polynomial.evaluate(p_eval);
            let d_eval = if i <= P_MAX_DEGREE {
                0.into()
            } else {
                let z_eval = PolynomialRG::interpolate_and_evaluate_zpoly(
                    0..P_MAX_DEGREE,
                    &i_as_field_element,
                );
                cp_eval.mul(
                    &z_eval
                        .inv()
                        .ok_or(anyhow!("unable to find inverse of Z evaluation"))
                        .unwrap(),
                )
            };

            original_evaluations.push(OriginalEvaluation::new(p_eval, d_eval));
        }

        let original_commitments_tree = MerkleTreeV2::new(20, &original_evaluations).unwrap();

        // We pseudo randomly generate a point for the column evaluation
        let x_c =
            pseudo_random_select_column_point(original_commitments_tree.root_hash(), TOTAL_POINTS);
        // We group the 4th root, i.e. all x's that have the same x^4 value, there are (p - 1)/4 + 1 = 262148 in total (over N - 1)
        let root_4_groups = group_roots_4(N - 1);
        // We evaluate g(x_c, x^4) over all the possible roots, it is a 250 degree polynomials
        let mut column_evaluations = vec![];
        for root_4_group in &root_4_groups {
            let y = root_4_group.four_exp();
            let evaluation = p.evaluate_as_binomial(x_c, y, 4);
            column_evaluations.push(ColumnEvaluation::new(evaluation));
        }
        let column_commitments_tree = MerkleTreeV2::new(19, &column_evaluations).unwrap();

        // Now we can generate the proof:
        // - original commitments / diagonal root hash,
        // - column commitments root hash,
        // - 300 random roots with their column evaluations and commitments in order to allow the verifier to low degree test the column at degree 250,
        // - 16 roots with their powers and their evaluations and commitments in order to generate 16 rows with 4 points each, roots must be chosen among the 300 random roots for the column in order to check consistence of the rows.

        let selected_column_points = pseudo_random_select_column_points(
            original_commitments_tree.root_hash(),
            &root_4_groups,
        );

        let mut original_commitments: Vec<FourRootOriginalCommitments<N>> = vec![];
        for (root_index, _) in selected_column_points.iter().take(16) {
            let root_group = &root_4_groups[*root_index as usize];
            let mut four_root_original_commitments: Vec<OriginalCommitment<N>> = vec![];
            for root in root_group.iter() {
                let selected_value = original_commitments_tree
                    .get(&point_to_selector(*root, 20))
                    .unwrap();
                four_root_original_commitments.push(OriginalCommitment::<N> {
                    point: *root,
                    evaluation: selected_value.value.clone(),
                    proof: selected_value.proof,
                });
            }
            original_commitments.push(four_root_original_commitments.try_into().unwrap());
        }

        let mut column_commitments: Vec<ColumnCommitment<N>> = vec![];
        for (root_index, root) in selected_column_points {
            let point_index_as_fe = FieldElementRG::from(root_index);
            let selected = column_commitments_tree
                .get(&point_to_selector(point_index_as_fe, 19))
                .unwrap();

            column_commitments.push(ColumnCommitment::<N> {
                point: root.exp(4),
                evaluation: selected.value.clone(),
                proof: selected.proof,
            });
        }

        let stark_proof = StarkProof {
            column_root: *column_commitments_tree.root_hash(),
            column: column_commitments,
            original_commitment_root: *original_commitments_tree.root_hash(),
            original: original_commitments,
        };

        // ######################
        // ###### Bob part ######
        // ######################

        let root_4_groups = group_roots_4(N - 1);
        let expected_x_c =
            pseudo_random_select_column_point(&stark_proof.original_commitment_root, TOTAL_POINTS);
        let expected_column_points = pseudo_random_select_column_points(
            &stark_proof.original_commitment_root,
            &root_4_groups,
        );

        assert_eq!(expected_x_c, x_c);

        /*
         * For the diagonal:
         *  - the i-th point must match the expected one from the prf,
         *  - in each 4-root:
         *      - the modular exponentiation by 4 must give the same result,
         *      - the merkle proof must be valid with respect to diagonal commitments,
         *      - the value of cp and d must be consistent
         *  - the 4 values plus the associated column value must form a degree <4 polynomial
         */
        for (i, four_root_original_commitment) in stark_proof.original.iter().enumerate() {
            // The i-th proof must be done with the i-th column root
            assert_eq!(
                four_root_original_commitment[0].point,
                expected_column_points[i].1
            );

            let expected_exp_4 = four_root_original_commitment[0].point.exp(4);

            let mut row_points = vec![];
            let mut row_values = vec![];
            for original_commitment in four_root_original_commitment {
                // Each sub point must give the same result under modular exp of 4
                assert_eq!(original_commitment.point.exp(4), expected_exp_4);
                // Merkle proof must be valid with respect to diagonal commitments
                assert!(
                    verify_proof(
                        &stark_proof.original_commitment_root,
                        &original_commitment.evaluation,
                        &original_commitment.proof
                    )
                    .is_ok()
                );

                let cp_eval = constraint_polynomial.evaluate(original_commitment.evaluation.p);
                if original_commitment.point.inner() < P_MAX_DEGREE {
                    // Below max degree, cp must be 0
                    assert_eq!(cp_eval, FieldElementRG::from(0));
                } else {
                    // Else, cp must be equal to z * d
                    let z_eval = PolynomialRG::interpolate_and_evaluate_zpoly(
                        0..P_MAX_DEGREE,
                        &original_commitment.point,
                    );
                    assert_eq!(z_eval.mul(&original_commitment.evaluation.d), cp_eval);
                }
                row_points.push(original_commitment.point);
                row_values.push(original_commitment.evaluation.p);
            }

            let matching_column_proof = &stark_proof.column[i];
            assert_eq!(matching_column_proof.point, expected_exp_4);

            row_points.push(expected_x_c);
            row_values.push(matching_column_proof.evaluation.v);
            let interpolated_polynomial =
                PolynomialRG::interpolate_from_coordinates(row_points, row_values).unwrap();
            assert!(interpolated_polynomial.degree() < 4);
        }

        let mut interpolation_points = vec![];
        let mut interpolation_values = vec![];
        for (i, column_commitment) in stark_proof.column.iter().enumerate() {
            // The i-th proof must be done with the i-th column root
            assert_eq!(column_commitment.point, expected_column_points[i].1.exp(4));
            // Merkle proof must be valid with respect to column commitments
            assert!(
                verify_proof(
                    &stark_proof.column_root,
                    &column_commitment.evaluation,
                    &column_commitment.proof
                )
                .is_ok()
            );
            interpolation_points.push(column_commitment.point);
            interpolation_values.push(column_commitment.evaluation.v);
        }
        let column_interpolated_polynomial =
            Polynomial::interpolate_from_coordinates(interpolation_points, interpolation_values)
                .unwrap();
        assert!(column_interpolated_polynomial.degree() < 250);
    }
}
