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

#[derive(Clone, Debug)]
pub struct StarkLeaf2<const N: u32> {
    pub v: PrimeFieldElement<N>,
    rep: [u8; 4],
}

impl<const N: u32> StarkLeaf2<N> {
    pub fn new(v: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 4];
        rep[0..4].clone_from_slice(&v.inner().to_le_bytes());
        Self { v, rep }
    }
}

impl<const N: u32> AsRef<[u8]> for StarkLeaf2<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u32> Default for StarkLeaf2<N> {
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
        merkletree_v2::{MerkleTreeV2, ValueWithProof, verify_proof},
        stark::polynomial::Polynomial,
    };

    const N: u32 = 1_048_589; // Chosen because: (p - 1) / 4 = 262147 and p > 2**20 (depth of the tree)
    const MINUS_ONE_SQRT_0: u32 = 38993; // sqrt(-1) (mod N)
    const MINUS_ONE_SQRT_1: u32 = 1009596; // sqrt(-1) (mod N)
    const P_MAX_DEGREE: u32 = 1_000; // 0 <= P(x) <= 9 for 0 <= x <= 999 => 1_000 constraints. One can always use Lagrange to build a valid 1_000 degree polynomial
    const TOTAL_POINTS: u32 = 1_048_576; // We take this in order to completely fill a Merkle tree of depth 20 while being lower than N

    type FieldElementRG = PrimeFieldElement<N>;
    type PolynomialRG = Polynomial<N>;

    fn generate_commitments_tree(
        p: &PolynomialRG,
        constraint_polynomial: &PolynomialRG,
        constrained_points: &[FieldElementRG],
    ) -> Result<MerkleTreeV2<StarkLeaf<N>>, anyhow::Error> {
        let mut values: Vec<StarkLeaf<N>> = Vec::with_capacity(TOTAL_POINTS as usize);
        for i in 0u32..TOTAL_POINTS {
            let i_as_field_element = FieldElementRG::from(i);

            let p_eval = p.evaluate(i_as_field_element);
            let cp_eval = constraint_polynomial.evaluate(p_eval);
            let d_eval = if i <= P_MAX_DEGREE {
                0.into()
            } else {
                let z_eval = PolynomialRG::interpolate_and_evaluate_from_roots_slice(
                    constrained_points,
                    &i_as_field_element,
                );
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

    fn derive_proofs<'a>(
        tree: &'a MerkleTreeV2<StarkLeaf<N>>,
        points: Vec<PrimeFieldElement<N>>,
    ) -> Result<Vec<(PrimeFieldElement<N>, ValueWithProof<'a, StarkLeaf<N>>)>, anyhow::Error> {
        let mut proofs: Vec<(PrimeFieldElement<N>, ValueWithProof<StarkLeaf<N>>)> = vec![];
        for point in points {
            let selector = tree.get(&point_to_selector(point, 20))?;
            proofs.push((point, selector))
        }
        Ok(proofs)
    }

    fn pseudo_random_select_points(root: &[u8; 32], modulus: u32) -> Vec<PrimeFieldElement<N>> {
        let mut points = Vec::with_capacity(16);
        for chunk in root.chunks(2) {
            let digest = Sha3_256::digest(chunk);
            let mut le_bytes = [0u8; 4];
            le_bytes.clone_from_slice(&digest[0..4]);
            let mut point_as_u32 = u32::from_le_bytes(le_bytes) % modulus;
            if point_as_u32 == 0 {
                point_as_u32 = 1;
            }
            points.push(PrimeFieldElement::<N>::from(point_as_u32));
        }
        points
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
        root_4_groups: &[[FieldElementRG; 4]],
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
                println!("selecting: {root_4_group_index}");
                if !selected_indices.contains(&root_4_group_index) {
                    let root_4_group = root_4_groups[root_4_group_index as usize];
                    if root_4_group.iter().all(|el| el.inner() < TOTAL_POINTS) {
                        points.push((root_4_group_index, root_4_group[0]));
                        selected_indices.insert(root_4_group_index);
                    } else {
                        println!("Invalid values: {root_4_group:?}\n");
                    }
                } else {
                    println!("Already selected\n")
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

    fn group_roots_4(max: u32) -> Vec<[FieldElementRG; 4]> {
        let zero = PrimeFieldElement::from(0);
        let mut elements = vec![[zero; 4]];

        let mut covered_elements = HashSet::new();

        let minus_one_sqrt_0 = FieldElementRG::from(MINUS_ONE_SQRT_0);
        let minus_one_sqrt_1 = FieldElementRG::from(MINUS_ONE_SQRT_1);

        for i in 1..max {
            let i_as_field_element = FieldElementRG::from(i);
            if !covered_elements.contains(&i_as_field_element) {
                let root_1 = i_as_field_element.mul(&minus_one_sqrt_0);
                let root_2 = i_as_field_element.mul(&minus_one_sqrt_1);
                let root_3 = i_as_field_element.neg();
                let roots = [i_as_field_element, root_1, root_2, root_3];
                covered_elements.insert(root_1);
                covered_elements.insert(root_2);
                covered_elements.insert(root_3);
                elements.push(roots);
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
            let expected = root_group[0].exp(4);
            for v in root_group.iter().skip(1).map(|el| el.exp(4)) {
                assert_eq!(v, expected);
            }
        }
    }

    struct StarkProof {
        pub column_root: [u8; 32],
        // Array of (point, evaluation, proof)
        pub column: Vec<(FieldElementRG, FieldElementRG, Vec<(bool, [u8; 32])>)>,
        pub diagonal_root: [u8; 32],
        // Arry of 4-tuple ([root, cp_root_eval d_root_eval, proof]) with the 4 4-roots
        pub diagonal: Vec<
            [(
                FieldElementRG,
                FieldElementRG,
                FieldElementRG,
                Vec<(bool, [u8; 32])>,
            ); 4],
        >,
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        // ###################
        // ###### Setup ######
        // ###################

        let constrained_points: Vec<FieldElementRG> =
            (0..P_MAX_DEGREE).map(FieldElementRG::from).collect(); // x(x - 1)...(x - P_MAX_DEGREE)

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

        // We generate commitments for the polynomial over 0..TOTAL_POINTS
        let commitments_tree =
            generate_commitments_tree(&p, &constraint_polynomial, &constrained_points).unwrap();

        // We pseudo randomly generate a point for the column evaluation
        let x_c = pseudo_random_select_column_point(commitments_tree.root_hash(), TOTAL_POINTS);
        // We group the 4th root, i.e. all x's that have the same x^4 value, there are (p - 1)/4 + 1 = 262148 in total (over N - 1)
        let root_4_groups = group_roots_4(N - 1);
        // We evaluate g(x_c, x^4) over all the possible roots, it is a 250 degree polynomials
        let mut column_evaluations = vec![];
        for root_4_group in &root_4_groups {
            let y = root_4_group[0].exp(4);
            let evaluation = p.evaluate_as_binomial(x_c, y, 4);
            column_evaluations.push(StarkLeaf2::new(evaluation));
        }
        let column_commitment_tree = MerkleTreeV2::new(19, &column_evaluations).unwrap();

        // Now we can generate the proof:
        // - diagonal commitments root hash,
        // - column commitments root hash,
        // - 300 random roots with their column evaluations and commitments in order to allow the verifier to low degree test the column at degree 250,
        // - 16 roots with their powers and their evaluations and commitments in order to generate 16 rows with 4 points each, roots must be chosen among the 300 random roots for the column in order to check consistence of the rows.

        let column_poins =
            pseudo_random_select_column_points(commitments_tree.root_hash(), &root_4_groups);

        let mut diagonal_proofs: Vec<
            [(
                FieldElementRG,
                FieldElementRG,
                FieldElementRG,
                Vec<(bool, [u8; 32])>,
            ); 4],
        > = vec![];
        for (root_index, _) in column_poins.iter().take(16) {
            let root_group = &root_4_groups[*root_index as usize];
            let mut group_proofs = vec![];
            for root in root_group {
                let selected_value = commitments_tree.get(&point_to_selector(*root, 20)).unwrap();
                let formatted_proof: Vec<(bool, [u8; 32])> = selected_value
                    .proof
                    .iter()
                    .map(|el| (el.0, *el.1))
                    .collect();
                group_proofs.push((
                    *root,
                    selected_value.value.cp,
                    selected_value.value.d,
                    formatted_proof,
                ));
            }
            diagonal_proofs.push(group_proofs.try_into().unwrap());
        }

        let mut column_proofs: Vec<(FieldElementRG, FieldElementRG, Vec<(bool, [u8; 32])>)> =
            vec![];
        for (root_index, root) in column_poins {
            let point_index_as_fe = FieldElementRG::from(root_index);
            let leaf = column_commitment_tree
                .get(&point_to_selector(point_index_as_fe, 19))
                .unwrap();

            let formatted_proof = leaf.proof.iter().map(|el| (el.0, *el.1)).collect();
            column_proofs.push((root, leaf.value.v, formatted_proof));
        }

        let _stark_proof = StarkProof {
            column_root: *column_commitment_tree.root_hash(),
            column: column_proofs,
            diagonal_root: *commitments_tree.root_hash(),
            diagonal: diagonal_proofs,
        };

        // DEPRECATED BELOW

        let selected_points =
            pseudo_random_select_points(commitments_tree.root_hash(), TOTAL_POINTS);
        let proofs = derive_proofs(&commitments_tree, selected_points).unwrap();

        // ######################
        // ###### Bob part ######
        // ######################

        let expected_points =
            pseudo_random_select_points(commitments_tree.root_hash(), TOTAL_POINTS);

        for ((point, value_with_proof), expected_point) in proofs.iter().zip(&expected_points) {
            assert_eq!(point, expected_point);
            if point.inner() < P_MAX_DEGREE {
                assert_eq!(value_with_proof.value.cp, 0.into());
            } else {
                let z_eval = PolynomialRG::interpolate_and_evaluate_from_roots_slice(
                    &constrained_points,
                    point,
                );
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
