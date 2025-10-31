use std::collections::HashSet;

use anyhow::anyhow;

use crate::{
    merkletree_v2::{MerkleProof, MerkleTreeV2, verify_proof},
    primefield::PrimeFieldElement,
    stark::{
        commitment::{Evaluation, commit},
        fri::{LowDegreeProof, generates_low_degree_proof, verify_low_degree_proof},
        polynomial::Polynomial,
        prf::pseudo_random_select_units_indices,
    },
};

const SPOT_CHECKS_COUNT: usize = 40;

/// Generates a stark proof for the following range check problem:
/// ```txt
/// The prover knows a polynomial `P` such that `0 <= P(x) <= 9 for 1 <= x <= MAX_DEGREE and P is of degree less than MAX_DEGREE`
/// ```
///
/// Let us define the additional properties of the system:
///     - N: the prime defining the prime field,
///     - C(x): the constraint polynomial defined as `C(x) = 0 for 0 <= x <= 9`, we take `C(x) = x(x - 1)...(x - 9)`,
///     - g: the generator of the prime field, i.e. any unit (all elements except 0) can be obtained as `a = g^k`.
///         The codebase will first organize the space as successive powers of this generator, hence we often refer to `unit index` (0 index is `g`, index 1 is `g^2`, etc...) instead of the unit directly,
///     - Z(x): the polynomial defined as `Z(x) = (x - 1)(x - 2)...(x - MAX_DEGREE)`.
///
/// With these definitions, our problem is now defined under the equation:
/// ```txt
/// C(P(x)) = Z(x) * D(x) (#1)
/// ```
/// Where the polynomial `D` has been introduced as the quotient of `C(P(x))` by `Z(x)`.
///
/// The proof generation is as follows:
///     1. we evaluate over the entire prime field the polynomials P and D.
///         Evaluations of D are made directly using the evaluations of `C(P)` and `Z`, we do not proceed by full interpolation.
///         Since `C(P)` is zero within the constrained interval, we don't have the associated `D` values and we put an arbitrary value.
///     2. evaluations of `D` and `P` are committed in two Merkle trees,
///     3. we select a number of spot checks. A spot check will be a query in the commitments of `P` and `D`.
///        The unit index, evaluation of P along with its Merkle proof, evaluation of D along with its Merkle proof define a spot check.
///        The indices are pseudo randomly chosen based on the root hash of the P commitments.
///     4. we generate the low degree proof that `P` is of degree less than `MAX_DEGREE`,
///     5. from equation (#1), the degree of `D` is at most `9 * MAX_DEGREE`, we generate the associated low degree proof.
///        The proof involves a lot of pseudo random selection of unit indices, in this proof, we need to make sure we don't select among the corrupted indices where evaluations of `D` do not make sense.
///
/// # Arguments
/// * `p` - the polynomial for which the proof is generated,
/// * `generator` - the generator of the prime field,
/// * `max_degree` - the maximum degree for the input polynomial
pub fn generate_stark_proof<const N: u64>(
    p: Polynomial<N>,
    generator: PrimeFieldElement<N>,
    max_degree: u64,
) -> Result<StarkProof<N>, anyhow::Error> {
    let units = derive_units(generator);

    let p_evals = p.fft_evaluate(&units);

    let mut p_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());
    let mut d_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());
    let mut invalid_d_evaluations_indices = HashSet::new();

    for (i, unit) in units.iter().enumerate() {
        let cp_eval = Polynomial::interpolate_and_evaluate_zpoly(0..10, &p_evals[i]);
        let d_eval = if unit.inner() <= max_degree {
            invalid_d_evaluations_indices.insert(i);
            // This value will not be checked so we can put any value we want here
            0.into()
        } else {
            let z_eval = Polynomial::interpolate_and_evaluate_zpoly(1..=max_degree, unit);
            cp_eval.mul(
                &z_eval
                    .inv()
                    .ok_or(anyhow!("unable to find inverse of Z evaluation"))?,
            )
        };

        p_evaluations.push(Evaluation::new(p_evals[i]));
        d_evaluations.push(Evaluation::new(d_eval));
    }

    let p_commitments_tree = commit(&p_evaluations)?;
    let d_commitments_tree = commit(&d_evaluations)?;

    let spot_checks = select_spot_checks(&p_commitments_tree, &d_commitments_tree, &units)
        .map_err(|e| e.context("selection of spot checks"))?;

    let p_low_degree_proof =
        generates_low_degree_proof(p_evaluations, p_commitments_tree, &units, max_degree, None)
            .map_err(|e| e.context("generation of low degree proof for p polynomial"))?;
    let d_low_degree_proof = generates_low_degree_proof(
        d_evaluations,
        d_commitments_tree,
        &units,
        9 * max_degree,
        Some(invalid_d_evaluations_indices),
    )
    .map_err(|e| e.context("generation of low degree proof for d polynomial"))?;

    Ok(StarkProof {
        spot_checks,
        p_low_degree_proof,
        d_low_degree_proof,
    })
}

fn select_spot_checks<const N: u64>(
    p_commitments_tree: &MerkleTreeV2<Evaluation<N>>,
    d_commitments_tree: &MerkleTreeV2<Evaluation<N>>,
    units: &[PrimeFieldElement<N>],
) -> Result<Vec<SpotCheck<N>>, anyhow::Error> {
    let spot_check_indices = pseudo_random_select_units_indices(
        p_commitments_tree.root_hash(),
        SPOT_CHECKS_COUNT,
        units.len(),
        &HashSet::new(),
    )?;
    let mut spot_checks = vec![];
    for unit_index in spot_check_indices {
        let p_commit = p_commitments_tree.select_commitment(unit_index)?;
        let d_commit = d_commitments_tree.select_commitment(unit_index)?;
        spot_checks.push(SpotCheck {
            unit_index,
            p_evaluation: p_commit.evaluation,
            p_proof: p_commit.proof,
            d_evaluation: d_commit.evaluation,
            d_proof: d_commit.proof,
        })
    }
    Ok(spot_checks)
}

/// Generates a stark proof for the following range check problem:
/// ```txt
/// The prover knows a polynomial `P` such that `0 <= P(x) <= 9 for 1 <= x <= MAX_DEGREE and P is of degree less than MAX_DEGREE`
/// ```
///
/// Let us define the additional properties of the system:
///     - N: the prime defining the prime field,
///     - C(x): the constraint polynomial defined as `C(x) = 0 for 0 <= x <= 9`, we take `C(x) = x(x - 1)...(x - 9)`,
///     - g: the generator of the prime field, i.e. any unit (all elements except 0) can be obtained as `a = g^k`.
///         The codebase will first organize the space as successive powers of this generator, hence we often refer to `unit index` (0 index is `g`, index 1 is `g^2`, etc...) instead of the unit directly,
///     - Z(x): the polynomial defined as `Z(x) = (x - 1)(x - 2)...(x - MAX_DEGREE)`.
///
/// With these definitions, our problem is now defined under the equation:
/// ```txt
/// C(P(x)) = Z(x) * D(x) (#1)
/// ```
/// Where the polynomial `D` has been introduced as the quotient of `C(P(x))` by `Z(x)`.
///
/// The verification is as follows:
///     - the low degree proof for `P` will be verified, it convinces the verifier, with good probability, that `P` is of degree less than `MAX_DEGREE`,
///     - the low degree proof for `D` will be verified, it convinces the verifier, with good probability, that `D` is of degree less than `9 * MAX_DEGREE`,
///     - the spot checks are used in verified in order to check that the relation #1 holds at the selected points.
/// The three checks together convince the verifier that the prover has a valid polynomial `P`.
///
/// # Arguments
/// * `stark_proof` - the stark proof to verify,
/// * `generator` - the generator of the prime field,
/// * `max_degree` - the maximum degree for the problem polynomial.
pub fn verify_stark_proof<const N: u64>(
    stark_proof: StarkProof<N>,
    generator: PrimeFieldElement<N>,
    max_degree: u64,
) -> Result<(), anyhow::Error> {
    let units = derive_units(generator);
    let mut invalid_d_evaluations_indices = HashSet::new();
    for (i, unit) in units.iter().enumerate() {
        if unit.inner() <= max_degree {
            invalid_d_evaluations_indices.insert(i);
        }
    }

    verify_low_degree_proof(&stark_proof.p_low_degree_proof, &units, max_degree, None)
        .map_err(|e| e.context("low degree test of p polynomial"))?;

    verify_low_degree_proof(
        &stark_proof.d_low_degree_proof,
        &units,
        9 * max_degree,
        Some(invalid_d_evaluations_indices),
    )
    .map_err(|e| e.context("low degree test of d polynomial"))?;

    verify_spot_checks(
        &stark_proof.p_low_degree_proof.original_commitments_root,
        &stark_proof.d_low_degree_proof.original_commitments_root,
        &stark_proof.spot_checks,
        max_degree,
        &units,
    )
    .map_err(|e| e.context("spot checks verification"))?;

    Ok(())
}

fn verify_spot_checks<const N: u64>(
    original_p_commitments_root: &[u8; 32],
    original_d_commitments_root: &[u8; 32],
    spot_checks: &[SpotCheck<N>],
    max_degree: u64,
    units: &[PrimeFieldElement<N>],
) -> Result<(), anyhow::Error> {
    let expected_spot_check_indices = pseudo_random_select_units_indices(
        original_p_commitments_root,
        SPOT_CHECKS_COUNT,
        units.len(),
        &HashSet::new(),
    )?;

    if spot_checks.len() != SPOT_CHECKS_COUNT {
        return Err(anyhow!(
            "invalid number of spot checks, got {}, expected {SPOT_CHECKS_COUNT}",
            spot_checks.len()
        ));
    }

    for (spot_check, expected_unit_index) in spot_checks.iter().zip(expected_spot_check_indices) {
        if spot_check.unit_index != expected_unit_index {
            return Err(anyhow!(
                "invalid unit index for spot check, got {}, expected {expected_unit_index}",
                spot_check.unit_index
            ));
        }

        verify_proof(
            original_p_commitments_root,
            &spot_check.p_evaluation,
            &spot_check.p_proof,
        )?;
        verify_proof(
            original_d_commitments_root,
            &spot_check.d_evaluation,
            &spot_check.d_proof,
        )?;

        let unit = &units[expected_unit_index];
        let cp_eval = Polynomial::interpolate_and_evaluate_zpoly(0..10, &spot_check.p_evaluation.v);
        if unit.inner() <= max_degree {
            if cp_eval != 0.into() {
                return Err(anyhow!(
                    "evaluation of constraint polynomial is not zero within constraint limits, got {cp_eval}"
                ));
            }
        } else {
            let z_eval = Polynomial::<N>::interpolate_and_evaluate_zpoly(1..=max_degree, unit);
            let rhs = z_eval.mul(&spot_check.d_evaluation.v);
            if rhs != cp_eval {
                return Err(anyhow!(
                    "original commitments do not respect the equation cp(p) = z * d, rhs is {rhs} and lhs is {cp_eval}"
                ));
            }
        }
    }

    Ok(())
}

pub struct StarkProof<const N: u64> {
    pub spot_checks: Vec<SpotCheck<N>>,
    pub p_low_degree_proof: LowDegreeProof<N>,
    pub d_low_degree_proof: LowDegreeProof<N>,
}

pub struct SpotCheck<const N: u64> {
    pub unit_index: usize,
    pub p_evaluation: Evaluation<N>,
    pub p_proof: MerkleProof,
    pub d_evaluation: Evaluation<N>,
    pub d_proof: MerkleProof,
}

fn derive_units<const N: u64>(generator: PrimeFieldElement<N>) -> Vec<PrimeFieldElement<N>> {
    let one = PrimeFieldElement::<N>::from(1);
    if generator == one {
        return vec![one];
    }
    let mut units = vec![one, generator];

    let mut power = generator.mul(&generator);
    while power != one {
        units.push(power);
        power = power.mul(&generator);
    }

    units
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::stark::polynomial::Polynomial;

    const N: u64 = 65_537; // One of the prime Fermat number, i.e. 2^(2^k) + 1
    const GENERATOR: u64 = 3;
    const P_MAX_DEGREE: u64 = 256; // 0 <= P(x) <= 9 for 1 <= x <= P_MAX_DEGREE

    #[test]
    fn test_power_of_4() {
        let root = PrimeFieldElement::<N>::from(GENERATOR);
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
        while a.is_multiple_of(4) {
            reduction_counter += 1;
            a /= 4;
        }

        assert!(reduction_counter > 4);

        let mut generator = PrimeFieldElement::<N>::from(GENERATOR);
        let mut old_len = derive_units(generator).len();

        for _ in 0..reduction_counter {
            generator = generator.exp(4);
            let new_len = derive_units(generator).len();
            assert_eq!(new_len, old_len / 4);
            old_len = new_len;
        }
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof_degree_0() {
        let p = Polynomial::<N>::new(vec![4.into()]);
        assert_eq!(p.degree(), 0);

        let stark_proof = generate_stark_proof(p, GENERATOR.into(), P_MAX_DEGREE).unwrap();
        verify_stark_proof(stark_proof, GENERATOR.into(), P_MAX_DEGREE).unwrap();
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof_degree_1023() {
        // We define a polynomial P with degree less than P_MAX_DEGREE and verifying 0 <= P(x) <= 9 for x between 1 and P_MAX_DEGREE
        // P will be 0 at 1 <= x < P_MAX_DEGREE, it will gives a (P_MAX_DEGREE - 1) polynomial that we're gonna scale down with the evaluation at x = P_MAX_DEGREE
        let mut p = Polynomial::<N>::interpolate_from_roots(
            (1..P_MAX_DEGREE)
                .map(PrimeFieldElement::<N>::from)
                .collect::<Vec<PrimeFieldElement<N>>>()
                .as_slice(),
        );
        p = p.mul_by_scalar(&p.evaluate(&P_MAX_DEGREE.into()).inv().unwrap());
        assert_eq!(p.degree(), P_MAX_DEGREE as usize - 1);

        let stark_proof = generate_stark_proof(p, GENERATOR.into(), P_MAX_DEGREE).unwrap();
        verify_stark_proof(stark_proof, GENERATOR.into(), P_MAX_DEGREE).unwrap();
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof_degree_1024() {
        let p = Polynomial::<N>::interpolate_from_roots(
            (1..(P_MAX_DEGREE + 1))
                .map(PrimeFieldElement::<N>::from)
                .collect::<Vec<PrimeFieldElement<N>>>()
                .as_slice(),
        );
        assert_eq!(p.degree(), P_MAX_DEGREE as usize);

        let stark_proof = generate_stark_proof(p, GENERATOR.into(), P_MAX_DEGREE).unwrap();
        assert!(verify_stark_proof(stark_proof, GENERATOR.into(), P_MAX_DEGREE).is_err());
    }

    #[test]
    #[cfg_attr(not(feature = "stark"), ignore)]
    fn test_stark_proof_invalid_polynomial() {
        let mut p = Polynomial::<N>::interpolate_from_roots(
            &(1..P_MAX_DEGREE)
                .map(PrimeFieldElement::<N>::from)
                .collect::<Vec<PrimeFieldElement<N>>>(),
        );
        p = p.mul_by_scalar(
            &p.evaluate(&P_MAX_DEGREE.into())
                .inv()
                .unwrap()
                .mul(&PrimeFieldElement::<N>::from(10)),
        );
        assert_eq!(p.degree(), P_MAX_DEGREE as usize - 1);

        let stark_proof = generate_stark_proof(p, GENERATOR.into(), P_MAX_DEGREE).unwrap();
        assert!(verify_stark_proof(stark_proof, GENERATOR.into(), P_MAX_DEGREE).is_err());
    }
}
