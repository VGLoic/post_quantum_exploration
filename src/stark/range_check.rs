use std::collections::HashSet;

use anyhow::anyhow;

use crate::{
    primefield::PrimeFieldElement,
    stark::{
        commitment::Evaluation,
        fri::{LowDegreeProof, low_degree_proof, verify_low_degree_proof},
        polynomial::Polynomial,
        prf::pseudo_random_select_units_indices,
    },
};

pub fn stark_proof<const N: u32>(
    p: Polynomial<N>,
    root_of_unity: PrimeFieldElement<N>,
    max_degree: u32,
) -> Result<StarkProof<N>, anyhow::Error> {
    let units = derive_units(root_of_unity);

    let mut p_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());
    let mut d_evaluations: Vec<Evaluation<N>> = Vec::with_capacity(units.len());

    let mut all_points = vec![];
    let mut all_values = vec![];
    let mut points = vec![];
    let mut values = vec![];

    let mut excluded_indices = HashSet::new();

    for (i, unit) in units.iter().enumerate() {
        let p_eval = p.evaluate(unit);
        let cp_eval = Polynomial::interpolate_and_evaluate_zpoly(0..10, &p_eval);
        let d_eval = if unit.inner() <= max_degree {
            excluded_indices.insert(i);
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

        if unit.inner() > max_degree {
            points.push(*unit);
            values.push(d_eval);
        }
        all_points.push(*unit);
        all_values.push(d_eval);

        p_evaluations.push(Evaluation::new(p_eval));
        d_evaluations.push(Evaluation::new(d_eval));
    }

    let p_low_degree_proof =
        low_degree_proof(p_evaluations, &units, max_degree, None, &excluded_indices)
            .map_err(|e| e.context("generation of low degree proof for p polynomial"))?;
    let d_low_degree_proof = low_degree_proof(
        d_evaluations,
        &units,
        9 * max_degree,
        Some(p_low_degree_proof.original_commitments_root),
        &excluded_indices,
    )
    .map_err(|e| e.context("generation of low degree proof for d polynomial"))?;

    Ok(StarkProof {
        p_low_degree_proof,
        d_low_degree_proof,
    })
}

pub fn verify_stark_proof<const N: u32>(
    stark_proof: StarkProof<N>,
    root_of_unity: PrimeFieldElement<N>,
    max_degree: u32,
) -> Result<(), anyhow::Error> {
    let units = derive_units(root_of_unity);
    let mut excluded_indices = HashSet::new();
    for (i, unit) in units.iter().enumerate() {
        if unit.inner() <= max_degree {
            excluded_indices.insert(i);
        }
    }

    verify_low_degree_proof(
        &stark_proof.p_low_degree_proof,
        &units,
        max_degree,
        None,
        &excluded_indices,
    )
    .map_err(|e| e.context("low degree test of p polynomial"))?;
    verify_low_degree_proof(
        &stark_proof.d_low_degree_proof,
        &units,
        9 * max_degree,
        Some(stark_proof.p_low_degree_proof.original_commitments_root),
        &excluded_indices,
    )
    .map_err(|e| e.context("low degree test of d polynomial"))?;

    let p_original_commitments = stark_proof
        .p_low_degree_proof
        .indirect_commitments
        .first()
        .map(|c| &c.diagonal_commitments)
        .ok_or(anyhow!("expected at least one indirect commitment"))?;
    let d_original_commitments = stark_proof
        .d_low_degree_proof
        .indirect_commitments
        .first()
        .map(|c| &c.diagonal_commitments)
        .ok_or(anyhow!("expected at least one indirect commitment"))?;

    let x_c_index = *pseudo_random_select_units_indices(
        &stark_proof.p_low_degree_proof.original_commitments_root,
        1,
        units.len(),
        &excluded_indices,
        None,
    )
    .first()
    .unwrap();
    let mut column_excluded_indices = HashSet::new();
    for i in &excluded_indices {
        column_excluded_indices.insert(i % (units.len() / 4));
    }
    let expected_row_unit_indices = pseudo_random_select_units_indices(
        &stark_proof.p_low_degree_proof.original_commitments_root,
        40,
        units.len() / 4,
        &column_excluded_indices,
        Some(x_c_index),
    );

    for ((p_row_commitments, d_row_commitments), expected_row_index) in p_original_commitments
        .iter()
        .zip(d_original_commitments)
        .zip(expected_row_unit_indices)
    {
        for (p_commitment, d_commitment) in p_row_commitments.iter().zip(d_row_commitments) {
            if p_commitment.unit_index != d_commitment.unit_index {
                return Err(anyhow!(
                    "non matching index for original data check, got index {} for p and index {} for d",
                    p_commitment.unit_index,
                    d_commitment.unit_index
                ));
            }
            let unit = &units[p_commitment.unit_index];
            let cp_eval =
                Polynomial::interpolate_and_evaluate_zpoly(0..10, &p_commitment.evaluation.v);
            if unit.inner() <= max_degree {
                if cp_eval != 0.into() {
                    return Err(anyhow!(
                        "evaluation of constraint polynomial is not zero within constraint limits, got {cp_eval}"
                    ));
                }
            } else {
                let z_eval = Polynomial::<N>::interpolate_and_evaluate_zpoly(1..=max_degree, unit);
                let rhs = z_eval.mul(&d_commitment.evaluation.v);
                if rhs != cp_eval {
                    return Err(anyhow!(
                        "original commitments do not respect the equation cp(p) = z * d, rhs is {rhs} and lfs is {cp_eval}"
                    ));
                }
            }
        }
        if p_row_commitments[0].unit_index != expected_row_index {
            return Err(anyhow!(
                "invalid value for unit index in row commitment, expected {expected_row_index}, got {}",
                p_row_commitments[0].unit_index
            ));
        }
    }

    Ok(())
}

pub struct StarkProof<const N: u32> {
    pub p_low_degree_proof: LowDegreeProof<N>,
    pub d_low_degree_proof: LowDegreeProof<N>,
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

    const N: u32 = 40_961; // Chosen because: (p - 1) / 4 = 268288 and larger than 1_000_000
    const ROOT_OF_UNITY: u32 = 3;
    // const N: u32 = 12_289; // Chosen because: (p - 1) / 4 = 268288 and larger than 1_000_000
    // const ROOT_OF_UNITY: u32 = 11;
    const P_MAX_DEGREE: u32 = 256; // 0 <= P(x) <= 9 for 1 <= x <= P_MAX_DEGREE
    // const N: u32 = 1_073_153; // Chosen because: (p - 1) / 4 = 268288 and larger than 1_000_000
    // const ROOT_OF_UNITY: u32 = 3;
    // const P_MAX_DEGREE: u32 = 1_024; // 0 <= P(x) <= 9 for 1 <= x <= P_MAX_DEGREE
    // REMIND ME: consider distinguishing max degree and constraint ranges

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

    // REMIND ME NEED TO DO SOMETHING?
    // #[test]
    // fn test_indices_selection() {
    //     let root = PrimeFieldElement::<N>::from(ROOT_OF_UNITY);
    //     let units = derive_units(root);

    //     let seed: [u8; 32] = rand::random();

    //     let mut excluded_indices = HashSet::new();
    //     for (i, unit) in units.iter().enumerate() {
    //         if unit.inner() <= P_MAX_DEGREE {
    //             excluded_indices.insert(i);
    //         }
    //     }

    //     let mut reduced_units = units;
    //     for reduction_level in 0..3 {
    //         let mut valid_indices_counter = 0;
    //         for (i, _) in reduced_units.iter().enumerate() {
    //             let is_excluded_index = excluded_indices.contains(&i);
    //             if !is_excluded_index {
    //                 valid_indices_counter += 1;
    //             }
    //         }

    //         // pseudo_random_select_units_indices(
    //         //     &seed,
    //         //     40,
    //         //     reduced_units.len(),
    //         //     reduction_level,
    //         //     &excluded_indices,
    //         // );
    //         println!(
    //             "Reduction #{reduction_level}:\n   Valid indices count: {valid_indices_counter} (among {} units)\n   Difference: {}\n####",
    //             reduced_units.len(),
    //             reduced_units.len() - valid_indices_counter,
    //         );
    //         assert!(valid_indices_counter >= 40);

    //         let mut newly_excluded_indices: HashSet<usize> = HashSet::new();
    //         for i in &excluded_indices {
    //             newly_excluded_indices.insert(i % (reduced_units.len() / 4));
    //         }
    //         reduced_units = reduced_units.into_iter().skip(3).step_by(4).collect();
    //         excluded_indices = newly_excluded_indices;
    //     }

    //     panic!("out");
    // }

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
    fn test_stark_proof_degree_0() {
        let p = Polynomial::<N>::new(vec![4.into()]);
        assert_eq!(p.degree(), 0);

        let stark_proof = stark_proof(p, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
        verify_stark_proof(stark_proof, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
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

        let stark_proof = stark_proof(p, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
        verify_stark_proof(stark_proof, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
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

        let stark_proof = stark_proof(p, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
        assert!(verify_stark_proof(stark_proof, ROOT_OF_UNITY.into(), P_MAX_DEGREE).is_err());
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

        let stark_proof = stark_proof(p, ROOT_OF_UNITY.into(), P_MAX_DEGREE).unwrap();
        assert!(verify_stark_proof(stark_proof, ROOT_OF_UNITY.into(), P_MAX_DEGREE).is_err());
    }
}
