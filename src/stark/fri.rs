use std::collections::HashSet;

use anyhow::anyhow;

use crate::{
    merkletree_v2::{MerkleTreeV2, verify_proof},
    primefield::PrimeFieldElement,
    stark::{
        commitment::{Commitment, Evaluation, commit},
        polynomial::Polynomial,
        prf::pseudo_random_select_units_indices,
    },
};

const DEGREE_THRESHOLD: u32 = 64;
const DIRECT_COMMITMENTS_COUNT: usize = 80;
const ROW_COUNT: usize = 40;

pub fn low_degree_proof<const N: u32>(
    original_values: Vec<Evaluation<N>>,
    original_commitments_tree: MerkleTreeV2<Evaluation<N>>,
    units: &[PrimeFieldElement<N>],
    max_degree: u32,
    seed: Option<[u8; 32]>,
    excluded_indices: &HashSet<usize>,
) -> Result<LowDegreeProof<N>, anyhow::Error> {
    let indirect_steps_count = derive_indirect_steps_count(max_degree)?;

    let mut diag_evaluations = original_values;
    let mut diagonal_commitments_tree = original_commitments_tree;
    let mut row_dimension = units.len();
    let mut seed = seed.unwrap_or(*diagonal_commitments_tree.root_hash());

    let original_commitments_root = *diagonal_commitments_tree.root_hash();

    let mut row_excluded_indices = excluded_indices.clone();

    let mut indirect_commitments: Vec<IndirectCommitment<N>> = vec![];
    for reduction_level in 0..(indirect_steps_count as usize) {
        let x_c_index = *pseudo_random_select_units_indices(
            &seed,
            1,
            row_dimension,
            &row_excluded_indices,
            None,
        )
        .first()
        .unwrap();
        let mut column_excluded_indices = HashSet::new();
        for i in &row_excluded_indices {
            column_excluded_indices.insert(i % (row_dimension / 4));
        }
        let row_indices = pseudo_random_select_units_indices(
            &seed,
            ROW_COUNT,
            row_dimension / 4,
            &column_excluded_indices,
            Some(x_c_index),
        );
        let x_c = &units[scale_index(x_c_index, reduction_level)];

        let mut column_evaluations = vec![];
        for unit_index in 0..(row_dimension / 4) {
            let mut row_points = vec![];
            let mut row_values = vec![];
            for i in 0..4 {
                let index = (unit_index + i * row_dimension / 4) % row_dimension;
                row_points.push(units[scale_index(index, reduction_level)]);
                row_values.push(diag_evaluations[index].v);
            }
            let row_interpolated_p =
                Polynomial::<N>::interpolate_from_coordinates(&row_points, &row_values).ok_or(
                    anyhow!("unable to interpolate row polynomial while evaluating column"),
                )?;
            column_evaluations.push(Evaluation::new(row_interpolated_p.evaluate(x_c)));
        }
        let column_commitments_tree = commit(&column_evaluations)?;

        let mut indirect_diagonal_commitments = vec![];
        let mut indirect_column_commitments = vec![];

        for row_index in row_indices {
            let row_diagonal_commitments = (0..4)
                .map(|i| {
                    let index = (row_index + i * row_dimension / 4) % row_dimension;
                    diagonal_commitments_tree.select_commitment(index)
                })
                .collect::<Result<Vec<Commitment<N>>, anyhow::Error>>()?;

            indirect_diagonal_commitments.push(
                row_diagonal_commitments.try_into().map_err(|_| {
                    anyhow!("unable to map diagonal commitments to a 4-length array")
                })?,
            );
            let associated_column_commitment =
                column_commitments_tree.select_commitment(row_index)?;
            indirect_column_commitments.push(associated_column_commitment);
        }

        indirect_commitments.push(IndirectCommitment {
            diagonal_commitments: indirect_diagonal_commitments,
            column_root: *column_commitments_tree.root_hash(),
            column_commitments: indirect_column_commitments,
        });

        row_dimension /= 4;
        diagonal_commitments_tree = column_commitments_tree;
        seed = *diagonal_commitments_tree.root_hash();
        diag_evaluations = column_evaluations;
        row_excluded_indices = column_excluded_indices;
    }

    let final_indices = pseudo_random_select_units_indices(
        diagonal_commitments_tree.root_hash(),
        DIRECT_COMMITMENTS_COUNT,
        row_dimension,
        &row_excluded_indices,
        None,
    );
    let direct_commitments = final_indices
        .into_iter()
        .map(|unit_index| diagonal_commitments_tree.select_commitment(unit_index))
        .collect::<Result<Vec<Commitment<N>>, anyhow::Error>>()?;

    Ok(LowDegreeProof {
        original_commitments_root,
        indirect_commitments,
        direct_commitments,
    })
}

pub fn verify_low_degree_proof<const N: u32>(
    proof: &LowDegreeProof<N>,
    units: &[PrimeFieldElement<N>],
    max_degree: u32,
    seed: Option<[u8; 32]>,
    excluded_indices: &HashSet<usize>,
) -> Result<(), anyhow::Error> {
    let indirect_steps_count = derive_indirect_steps_count(max_degree)?;

    let mut seed: [u8; 32] = seed.unwrap_or(proof.original_commitments_root);
    let mut diag_root = proof.original_commitments_root;

    let mut row_dimension = units.len();
    let mut row_excluded_indices = excluded_indices.clone();

    if indirect_steps_count as usize != proof.indirect_commitments.len() {
        return Err(anyhow!(
            "invalid number of indirect steps, expected {indirect_steps_count}, got {}",
            proof.indirect_commitments.len()
        ));
    }

    for (reduction_level, indirect_commitment) in proof.indirect_commitments.iter().enumerate() {
        let x_c_index = *pseudo_random_select_units_indices(
            &seed,
            1,
            row_dimension,
            &row_excluded_indices,
            None,
        )
        .first()
        .unwrap();
        let mut column_excluded_indices = HashSet::new();
        for i in &row_excluded_indices {
            column_excluded_indices.insert(i % (row_dimension / 4));
        }
        let row_indices = pseudo_random_select_units_indices(
            &seed,
            ROW_COUNT,
            row_dimension / 4,
            &column_excluded_indices,
            Some(x_c_index),
        );
        let x_c = &units[scale_index(x_c_index, reduction_level)];

        if indirect_commitment.column_commitments.len()
            != indirect_commitment.diagonal_commitments.len()
        {
            return Err(anyhow!(
                "indirect commitment is invalid, column and diagonal commitments length do not match, got column: {}, diag: {}",
                indirect_commitment.column_commitments.len(),
                indirect_commitment.diagonal_commitments.len()
            ));
        }

        for ((diag_row_commitments, associated_column_commitment), expected_row_index) in
            indirect_commitment
                .diagonal_commitments
                .iter()
                .zip(&indirect_commitment.column_commitments)
                .zip(row_indices)
        {
            let expected_four_exp_value =
                &units[scale_index(expected_row_index, reduction_level)].exp(4);

            let row_unit =
                &units[scale_index(associated_column_commitment.unit_index, reduction_level + 1)];
            if row_unit != expected_four_exp_value {
                return Err(anyhow!(
                    "invalid unit of the column commitment, expected {expected_four_exp_value}, got {row_unit}"
                ));
            }

            let mut row_points = vec![];
            let mut row_values = vec![];
            for row_commitment in diag_row_commitments {
                verify_proof(
                    &diag_root,
                    &row_commitment.evaluation,
                    &row_commitment.proof,
                )?;

                let unit = &units[scale_index(row_commitment.unit_index, reduction_level)];
                if &unit.exp(4) != expected_four_exp_value {
                    return Err(anyhow!(
                        "the power of 4 of the commitment unit point must be equal to the column unit"
                    ));
                }

                row_points.push(*unit);
                row_values.push(row_commitment.evaluation.v);
            }

            verify_proof(
                &indirect_commitment.column_root,
                &associated_column_commitment.evaluation,
                &associated_column_commitment.proof,
            )?;

            row_points.push(*x_c);
            row_values.push(associated_column_commitment.evaluation.v);

            let interpolated_p =
                Polynomial::<N>::interpolate_from_coordinates(&row_points, &row_values)
                    .ok_or(anyhow!("failed to interpolate row polynomial"))?;

            if interpolated_p.degree() >= 4 {
                return Err(anyhow!(
                    "interpolated row polynomial has a degree larger or equal to 4, got: {}",
                    interpolated_p.degree()
                ));
            }
        }

        row_dimension /= 4;
        diag_root = indirect_commitment.column_root;
        seed = indirect_commitment.column_root;
        row_excluded_indices = column_excluded_indices;
    }

    let final_indices = pseudo_random_select_units_indices(
        &diag_root,
        DIRECT_COMMITMENTS_COUNT,
        row_dimension,
        &row_excluded_indices,
        None,
    );

    if final_indices.len() != proof.direct_commitments.len() {
        return Err(anyhow!(
            "invalid direct commitments, expected {DIRECT_COMMITMENTS_COUNT}, got {}",
            proof.direct_commitments.len()
        ));
    }
    let mut points = vec![];
    let mut values = vec![];

    for (direct_commitment, expected_unit_index) in
        proof.direct_commitments.iter().zip(final_indices)
    {
        verify_proof(
            &diag_root,
            &direct_commitment.evaluation,
            &direct_commitment.proof,
        )?;

        if expected_unit_index != direct_commitment.unit_index {
            return Err(anyhow!(
                "non matching unit index, expected {}, got {}",
                expected_unit_index,
                direct_commitment.unit_index
            ));
        }

        points.push(
            units[scale_index(
                direct_commitment.unit_index,
                proof.indirect_commitments.len(),
            )],
        );
        values.push(direct_commitment.evaluation.v);
    }

    let interpolated_p = Polynomial::<N>::interpolate_from_coordinates(&points, &values)
        .ok_or(anyhow!("failed to interpolate from direct commitments"))?;

    if interpolated_p.degree() >= DEGREE_THRESHOLD as usize {
        return Err(anyhow!(
            "direct commitments lead to an interpolation of degree higher than {DEGREE_THRESHOLD}, got: {}",
            interpolated_p.degree()
        ));
    }

    Ok(())
}

pub struct LowDegreeProof<const N: u32> {
    pub original_commitments_root: [u8; 32],
    pub indirect_commitments: Vec<IndirectCommitment<N>>,
    pub direct_commitments: Vec<Commitment<N>>,
}

pub struct IndirectCommitment<const N: u32> {
    // Vec of diagonal commitments, each element contains 4 values over a 4-root in order to form a row
    pub diagonal_commitments: Vec<[Commitment<N>; 4]>,
    pub column_root: [u8; 32],
    // Vec of column commitments, each element contains 4 values over a 4-root, only the first one should be used in order to match the row from the diagonal commitment above
    pub column_commitments: Vec<Commitment<N>>,
}

fn derive_indirect_steps_count(max_degree: u32) -> Result<u32, anyhow::Error> {
    if max_degree == 0 {
        return Ok(0);
    }
    let mut degree = max_degree;
    let mut i = 0;
    while degree > DEGREE_THRESHOLD {
        if degree % 4 != 0 {
            return Err(anyhow!(
                "max_degree must be a divisble by 4 until degree threshold {DEGREE_THRESHOLD} is reached"
            ));
        }
        degree /= 4;
        i += 1;
    }

    Ok(i)
}

fn scale_index(index: usize, reduction_level: usize) -> usize {
    let mut result = index;
    for _ in 0..reduction_level {
        result = 4 * result + 3;
    }
    result
}
