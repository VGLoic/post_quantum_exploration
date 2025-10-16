/*
 * Our goal is to implement a stark proof for the knowledge of a polynomial P of degree less than 1_000_000 such that 0 <= P(x) <= 9 for x between 1 and 1_000_000.
 */

use anyhow::anyhow;

use crate::{
    merkletree_v2::verify_proof,
    primefield::PrimeFieldElement,
    stark::{
        commitment::{Commitment, Evaluation, build_tree},
        fri::{LowDegreeProof, low_degree_proof, verify_low_degree_proof},
        polynomial::Polynomial,
        prf::pseudo_random_select_units_indices,
    },
};

const P_MAX_DEGREE: u32 = 1_000;

pub fn stark_proof<const N: u32>(
    p: Polynomial<N>,
    root_of_unity: PrimeFieldElement<N>,
) -> Result<StarkProof<N>, anyhow::Error> {
    let units = derive_units(root_of_unity);

    // We evaluate the polynomial over all the available units
    let mut p_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());
    let mut d_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());

    for unit in &units {
        let p_eval = p.evaluate(*unit);
        let cp_eval = Polynomial::interpolate_and_evaluate_zpoly(0..10, &p_eval);
        let d_eval = if unit.inner() <= P_MAX_DEGREE {
            0.into()
        } else {
            let z_eval = Polynomial::interpolate_and_evaluate_zpoly(1..(P_MAX_DEGREE + 1), unit);
            cp_eval.mul(
                &z_eval
                    .inv()
                    .ok_or(anyhow!("unable to find inverse of Z evaluation"))?,
            )
        };

        p_evaluations.push(Evaluation::new(p_eval));
        d_evaluations.push(Evaluation::new(d_eval));
    }

    let p_commitments_tree = build_tree(&p_evaluations)?;
    let d_commitments_tree = build_tree(&d_evaluations)?;

    let original_commitments_indices =
        pseudo_random_select_units_indices(p_commitments_tree.root_hash(), 40, units.len());
    let mut p_commitments = vec![];
    let mut d_commitments = vec![];
    for index in original_commitments_indices {
        p_commitments.push(p_commitments_tree.select_commitment(index)?);
        d_commitments.push(d_commitments_tree.select_commitment(index)?);
    }

    let original_commitment = OriginalCommitment {
        p_root: *p_commitments_tree.root_hash(),
        p_commitments,
        d_root: *d_commitments_tree.root_hash(),
        d_commitments,
    };

    let low_degree_proof = low_degree_proof(p, p_commitments_tree, units, P_MAX_DEGREE)?;

    Ok(StarkProof {
        low_degree_proof,
        original_commitment,
    })
}

pub fn verify_stark_proof<const N: u32>(
    stark_proof: StarkProof<N>,
    root_of_unity: PrimeFieldElement<N>,
) -> Result<(), anyhow::Error> {
    let units = derive_units(root_of_unity);

    let expected_number_of_original_commitments = 40;
    let expected_unit_indices = pseudo_random_select_units_indices(
        &stark_proof.original_commitment.p_root,
        expected_number_of_original_commitments,
        units.len(),
    );

    if stark_proof.original_commitment.d_commitments.len()
        != stark_proof.original_commitment.p_commitments.len()
    {
        return Err(anyhow!(
            "different number of commitments for `d` (#{}) and `p` (#{}) polynomials",
            stark_proof.original_commitment.d_commitments.len(),
            stark_proof.original_commitment.p_commitments.len()
        ));
    }

    if stark_proof.original_commitment.d_commitments.len() < expected_number_of_original_commitments
    {
        return Err(anyhow!(
            "not enough commitments for the original data, expected at least {expected_number_of_original_commitments}, got {}",
            stark_proof.original_commitment.d_commitments.len()
        ));
    }

    for ((p_commitment, d_commitment), expected_unit_index) in stark_proof
        .original_commitment
        .p_commitments
        .iter()
        .zip(&stark_proof.original_commitment.d_commitments)
        .zip(expected_unit_indices)
    {
        if p_commitment.unit_index != d_commitment.unit_index {
            return Err(anyhow!(
                "different unit used for evaluation of `p` (got index: {}) and `d` (got index: {}) polynomials",
                p_commitment.unit_index,
                d_commitment.unit_index
            ));
        }
        if p_commitment.unit_index != expected_unit_index {
            return Err(anyhow!(
                "evaluated unit is not the expected one, got index {}, expected index {expected_unit_index}",
                p_commitment.unit_index
            ));
        }

        verify_proof(
            &stark_proof.low_degree_proof.original_commitments_root,
            &p_commitment.evaluation,
            &p_commitment.proof,
        )?;
        verify_proof(
            &stark_proof.original_commitment.d_root,
            &d_commitment.evaluation,
            &d_commitment.proof,
        )?;

        let unit = &units[expected_unit_index];

        let cp_eval = Polynomial::interpolate_and_evaluate_zpoly(0..10, &p_commitment.evaluation.v);
        if unit.inner() <= P_MAX_DEGREE {
            if cp_eval != 0.into() {
                return Err(anyhow!(
                    "evaluation of constraint polynomial is not zero within constraint limits, got {cp_eval}"
                ));
            }
        } else {
            let z_eval =
                Polynomial::<N>::interpolate_and_evaluate_zpoly(1..(P_MAX_DEGREE + 1), unit);
            let rhs = z_eval.mul(&d_commitment.evaluation.v);
            if rhs != cp_eval {
                return Err(anyhow!(
                    "original commitments do not respect the equation cp(p) = z * d, rhs is {rhs} and lfs is {cp_eval}"
                ));
            }
        }
    }

    verify_low_degree_proof(stark_proof.low_degree_proof, units, P_MAX_DEGREE)?;

    Ok(())
}

pub struct StarkProof<const N: u32> {
    pub original_commitment: OriginalCommitment<N>,
    pub low_degree_proof: LowDegreeProof<N>,
}

pub struct OriginalCommitment<const N: u32> {
    pub p_root: [u8; 32],
    pub d_root: [u8; 32],
    pub p_commitments: Vec<Commitment<N>>,
    pub d_commitments: Vec<Commitment<N>>,
}

fn derive_units<const N: u32>(root: PrimeFieldElement<N>) -> Vec<PrimeFieldElement<N>> {
    let mut units = vec![root];

    let mut power = root.mul(&root);
    while power != root {
        units.push(power);
        power = power.mul(&root);
    }

    units
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::stark::polynomial::Polynomial;

    const N: u32 = 1_073_153; // Chosen because: (p - 1) / 4 = 268288 and larger than 1_000_000
    const ROOT_OF_UNITY: u32 = 3;
    const P_MAX_DEGREE: u32 = 1_000; // 0 <= P(x) <= 9 for 1 <= x <= 1_000 => 1_000 constraints. One can always use Lagrange to build a valid 1_000 degree polynomial // REMIND ME

    type FieldElementRG = PrimeFieldElement<N>;
    type PolynomialRG = Polynomial<N>;

    #[test]
    fn test_power_of_4() {
        let root = PrimeFieldElement::<N>::from(ROOT_OF_UNITY);
        let units = derive_units(root);
        for (i, unit) in units.iter().enumerate() {
            let expected_value = unit.exp(4);
            assert_eq!(
                units[(i + units.len() / 4) % units.len()].exp(4),
                expected_value
            );
            assert_eq!(
                units[(i + 2 * units.len() / 4) % units.len()].exp(4),
                expected_value
            );
            assert_eq!(
                units[(i + 3 * units.len() / 4) % units.len()].exp(4),
                expected_value
            );
        }
    }

    #[test]
    fn test_fri_friendly() {
        let mut reduction_counter = 0;
        let mut a = N - 1;
        while a % 4 == 0 {
            reduction_counter += 1;
            a /= 4;
        }

        assert!(reduction_counter > 4);

        let mut root = PrimeFieldElement::<N>::from(ROOT_OF_UNITY);
        let mut old_len = derive_units(root).len();

        for _ in 0..reduction_counter {
            root = root.exp(4);
            let new_len = derive_units(root).len();
            assert_eq!(new_len, old_len / 4);
            old_len = new_len;
        }
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof() {
        // We define a polynomial P with degree less than P_MAX_DEGREE and verifying 0 <= P(x) <= 9 for x between 1 and P_MAX_DEGREE
        // P will be 0 at 1 <= x < P_MAX_DEGREE, it will gives a (P_MAX_DEGREE - 1) polynomial that we're gonna scale down with the evaluation at x = P_MAX_DEGREE
        let mut p = PolynomialRG::interpolate_from_roots(
            (1..P_MAX_DEGREE).map(FieldElementRG::from).collect(),
        );
        p = p.mul_by_scalar(p.evaluate(P_MAX_DEGREE.into()).inv().unwrap());
        if p.degree() >= P_MAX_DEGREE as usize {
            panic!("Invalid degree, got {}", p.degree());
        }

        let stark_proof = stark_proof(p, ROOT_OF_UNITY.into()).unwrap();
        verify_stark_proof(stark_proof, ROOT_OF_UNITY.into()).unwrap();
    }
}
