use std::collections::HashSet;

use anyhow::anyhow;

use crate::{
    merkletree_v2::{MerkleTreeV2, verify_proof},
    primefield::PrimeFieldElement,
    stark::{
        commitment::{Commitment, Evaluation, commit},
        polynomial::Polynomial,
        prf::{pseudo_random_select_unit_index, pseudo_random_select_units_indices},
    },
};

const DEGREE_THRESHOLD: u64 = 64;
const DIRECT_COMMITMENTS_COUNT: usize = 80;
const ROW_COUNT: usize = 40;

/// Generates a low degree proof for a set of a polynomial evaluations
///
/// I.e. generates a proof that the set of evaluations is associated, with high probability, to a polynomial of degree less than an input.
///
/// This function uses a [FRI](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) like algorithm.
///
/// # Arguments
///
/// * `original_values` - the list of evaluations of the polynomial for which we provide the proof, the evaluations must be ordered as the successive powers of the generator of the prime field,
/// * `original_commitments_tree` - the Merkle tree formed from the list of evaluations, the coherence with the original evaluations input is not checked in the function,
/// * `units` - the list of units, ordered as successive powers of the generator of the prime field,
/// * `max_degree` - the degree under which the polynomial must belong,
/// * `excluded_indices` - a list of unit indices that should not be selected during proof generation.
///
/// # Implementation
///
/// For a more detailed version, please consult the README of the module.
///
/// Let us state two important remarks:
/// 1. this algorithm will only work with prime field containing a sufficiently large subgroup of order 2^k.
///    More concretely, for `N` the prime, we need prime field such that `N - 1` is divisible by the factor `4` as much as possible.
///    This is a required ingredient for the process of `space reduction` during the indirect proofs generation step.
/// 2. let us define `g` the generator of our prime field. The generator is used in order to derive the units as successive power of `g`.
///    I.e. `units = [g, g^2, g^3, ... , g^(N-1)]`.
///    The input evaluations are done over these "ordered" units.
///
/// The function starts by preparing for the recursive part:
/// - the number of indirect steps is pre-computed,
/// - the diagonal evaluations are instantiated as the original evaluations,
/// - the diagonal commitments are instantiated as the original commitments, the function does not check that these commitments match with the evaluations,
/// - the row dimension is instantiated as the number of units,
/// - the excluded indices along the row is instantiated with the provided list, if any,
/// - the original commitment root hash is saved for the output proof.
///
/// ## Indirect steps
///
/// At each step:
/// - the diagonal commitment root hash is used in order to pseudo randomly select the unit describing the column: `x_c`,
/// - the diagonal polynomial is evaluated as a binomial polynomial along the column `x = x_c`, for each unit in the column:
///     - the four diagonal evaluations on this row are selected using their indices: `column_unit_index + i * (N-1)/4, i = 0, 1, 2, 3`,
///     - the polynomial of the row is interpolated using Lagrange interpolation,
///     - the polynomial is evaluated at `x = x_c`,
/// - the column evaluations are committed,
/// - the excluded indices of the column are computed from the ones of the row: if a column index is on a row with an excluded index, it is excluded.
///   In practice, it means mapping `index -> index % (row_dimension / 4)`,
/// - the root of the column commit is used in order to pseudo randomly select indices of the rows that will be added to the proof for this step.
///   Indices must not be in the list of exluded ones. An additional index is taken in case we select the row crossing the diagonal at `x = x_c`, if this is the case, we remove this row and use the additional one instead,
/// - for each selected row, we select the four diagonal evaluations and their Merkle roots, we also select the associated column evaluation and its Merkle proof, these five commitments form the proof for the row of the indirect step,
/// - the next step is prepared:
///     - the row dimension becomes the column dimension, i.e. it is divided by 4,
///     - the diagonal evaluations and commitments are set to the column ones,
///     - the row excluded indices are set to the previously computed ones for the column
///
/// We end up with a list of indirect proof, each defined as a list of proof for a row, each defined as five proven evaluations:
/// - 4 coming from the diagonal evaluations,
/// - 1 coming from the column evaluation.
///
/// ## Final and direct step
///
/// - the final indices are pseudo randomly selected from the lastly defined diagonal commitments, these indices must not contain any lastly defined excluded indices,
/// - for each index, we select the evaluation with its Merkle proof.
///
/// ## Proof content
///
/// The proof is simply build as:
/// - the list of indirect proofs,
/// - the direct commitments,
/// - the root hash of the original commitments.
pub fn generates_low_degree_proof<const N: u64>(
    original_values: Vec<Evaluation<N>>,
    original_commitments_tree: MerkleTreeV2<Evaluation<N>>,
    units: &[PrimeFieldElement<N>],
    max_degree: u64,
    excluded_indices: Option<HashSet<usize>>,
) -> Result<LowDegreeProof<N>, anyhow::Error> {
    let indirect_steps_count = derive_indirect_steps_count(max_degree)?;

    let mut diag_evaluations = original_values;
    let mut diagonal_commitments_tree = original_commitments_tree;
    let mut row_dimension = units.len();

    let original_commitments_root = *diagonal_commitments_tree.root_hash();

    let mut row_excluded_indices = excluded_indices.unwrap_or_default();

    let mut indirect_commitments: Vec<IndirectCommitment<N>> = vec![];
    for reduction_level in 0..(indirect_steps_count as usize) {
        let x_c_index =
            pseudo_random_select_unit_index(diagonal_commitments_tree.root_hash(), row_dimension);

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
                Polynomial::<N>::interpolate_from_coordinates(&row_points, &row_values)
                    .map_err(|e| e.context("interpolation of row polynomial"))?;
            column_evaluations.push(Evaluation::new(row_interpolated_p.evaluate(x_c)));
        }
        let column_commitments_tree = commit(&column_evaluations)?;

        let mut column_excluded_indices = HashSet::new();
        for i in &row_excluded_indices {
            column_excluded_indices.insert(i % (row_dimension / 4));
        }
        let row_indices_with_additional_unit = pseudo_random_select_units_indices(
            column_commitments_tree.root_hash(),
            ROW_COUNT + 1,
            row_dimension / 4,
            &column_excluded_indices,
        )?;
        let mut indirect_diagonal_commitments = vec![];
        let mut indirect_column_commitments = vec![];
        let excluded_row_index = x_c_index % (row_dimension / 4);
        for row_index in row_indices_with_additional_unit.iter().take(ROW_COUNT) {
            let mut row_index = *row_index;
            // We don't want to select the unique row that cross with the selected column
            // otherwise we will end up with only four points on the line and not five
            if row_index == excluded_row_index {
                row_index = row_indices_with_additional_unit[ROW_COUNT];
            }

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
        diag_evaluations = column_evaluations;
        row_excluded_indices = column_excluded_indices;
    }

    let final_indices = pseudo_random_select_units_indices(
        diagonal_commitments_tree.root_hash(),
        DIRECT_COMMITMENTS_COUNT,
        row_dimension,
        &row_excluded_indices,
    )?;
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

/// Verifies a low degree proof
///
/// I.e. verify that a proof matches, with high probability, with a polynomial of degree less than an input.
///
/// This function uses a [FRI](https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html) like algorithm.
///
/// # Arguments
///
/// * `proof` - the low degree proof. It contains:
///     - the original commitments root hash,
///     - the list of indirect commitments,
///     - the list of direct commitments,
/// * `units` - the list of units, ordered as successive powers of the generator of the prime field,
/// * `max_degree` - the degree under which the polynomial must belong,
/// * `excluded_indices` - a list of unit indices that should not be selected during proof generation.
///
/// # Implementation
///
/// For a more detailed version, please consult the README of the module.
///
/// Let us state two important remarks:
/// 1. this algorithm will only work with prime field containing a sufficiently large subgroup of order 2^k.
///    More concretely, for `N` the prime, we need prime field such that `N - 1` is divisible by the factor `4` as much as possible.
///    This is a required ingredient for the process of `space reduction` during the indirect proofs generation step.
/// 2. let us define `g` the generator of our prime field. The generator is used in order to derive the units as successive power of `g`.
///    I.e. `units = [g, g^2, g^3, ... , g^(N-1)]`.
///    The input evaluations are done over these "ordered" units.
///
/// The function starts by preparing for the verification of the recursive/indirect part:
/// - the number of expected indirect steps is pre-computed,
/// - the diagonal commitment root hash is instantiated with the one of the original commitments,
/// - the row dimension is instantiated as the number of units,
/// - the excluded indices along the row is instantiated with the provided list, if any,
///
/// ## Verification of indirect steps
///
/// For each step:
/// - the diagonal commitment root hash is used in order to pseudo randomly select the unit describing the column: `x_c`,
/// - the excluded indices of the column are computed from the ones of the row: if a column index is on a row with an excluded index, it is excluded.
///   In practice, it means mapping `index -> index % (row_dimension / 4)`,
/// - the root of the column commit is used in order to pseudo randomly select indices of the rows that will be added to the proof for this step.
///   Indices must not be in the list of exluded ones. An additional index is taken in case we select the row crossing the diagonal at `x = x_c`, if this is the case, we remove this row and use the additional one instead,
/// - for each row:
///     - the units coming from the diagonal commits are verified to be on the same row as the expected one,
///     - the Merkle proofs of the diagonal commits are verified against the diagonal commitment root hash,
///     - the unit of the associated column commit is verified against the expected row,
///     - the Merkle proof of the associated column commit is verified against the column commitment root hash,
///     - the four evaluations of the diagonal at specified units added with the column evaluation at `x_c` are used for a polynomial interpolation,
///     - the interpolated polynomial must be of degree less than 4,
/// - the next step is prepared:
///     - the row dimension becomes the column dimension, i.e. it is divided by 4,
///     - the diagonal evaluations and commitments are set to the column ones,
///     - the row excluded indices are set to the previously computed ones for the column
///
/// ## Verification of the direct step
///
/// - the expected final indices are pseudo randomly selected from the lastly defined diagonal commitments, these indices must not contain any lastly defined excluded indices,
/// - for each index:
///     - the index is verified against the expected one,
///     - the Merkle proof is verified against the lastly defined diagonal commitment root hash,
/// - the list of evaluations is interpolated, the resulting polynomial must be of degree less than the set threshold.
pub fn verify_low_degree_proof<const N: u64>(
    proof: &LowDegreeProof<N>,
    units: &[PrimeFieldElement<N>],
    max_degree: u64,
    excluded_indices: Option<HashSet<usize>>,
) -> Result<(), anyhow::Error> {
    let indirect_steps_count = derive_indirect_steps_count(max_degree)?;

    let mut diag_root = proof.original_commitments_root;

    let mut row_dimension = units.len();
    let mut row_excluded_indices = excluded_indices.unwrap_or_default();

    if indirect_steps_count as usize != proof.indirect_commitments.len() {
        return Err(anyhow!(
            "invalid number of indirect steps, expected {indirect_steps_count}, got {}",
            proof.indirect_commitments.len()
        ));
    }

    for (reduction_level, indirect_commitment) in proof.indirect_commitments.iter().enumerate() {
        let x_c_index = pseudo_random_select_unit_index(&diag_root, row_dimension);
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

        let mut column_excluded_indices = HashSet::new();
        for i in &row_excluded_indices {
            column_excluded_indices.insert(i % (row_dimension / 4));
        }
        let row_indices_with_additional_unit = pseudo_random_select_units_indices(
            &indirect_commitment.column_root,
            ROW_COUNT + 1,
            row_dimension / 4,
            &column_excluded_indices,
        )?;
        let excluded_row_index = x_c_index % (row_dimension / 4);
        for ((diag_row_commitments, associated_column_commitment), expected_row_index) in
            indirect_commitment
                .diagonal_commitments
                .iter()
                .zip(&indirect_commitment.column_commitments)
                .zip(&row_indices_with_additional_unit)
        {
            let mut expected_row_index = *expected_row_index;
            // We don't want to select the unique row that cross with the selected column
            // otherwise we will end up with only four points on the line and not five
            if expected_row_index == excluded_row_index {
                expected_row_index = row_indices_with_additional_unit[ROW_COUNT];
            }
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
                    .map_err(|e| e.context("interpolation of row polynomial"))?;

            if interpolated_p.degree() >= 4 {
                return Err(anyhow!(
                    "interpolated row polynomial has a degree larger or equal to 4, got: {}",
                    interpolated_p.degree()
                ));
            }
        }

        row_dimension /= 4;
        diag_root = indirect_commitment.column_root;
        row_excluded_indices = column_excluded_indices;
    }

    let final_indices = pseudo_random_select_units_indices(
        &diag_root,
        DIRECT_COMMITMENTS_COUNT,
        row_dimension,
        &row_excluded_indices,
    )?;

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
        .map_err(|e| e.context("interpolation of polynomial for direct check"))?;

    if interpolated_p.degree() >= DEGREE_THRESHOLD as usize {
        return Err(anyhow!(
            "direct commitments lead to an interpolation of degree higher than {DEGREE_THRESHOLD}, got: {}",
            interpolated_p.degree()
        ));
    }

    Ok(())
}

pub struct LowDegreeProof<const N: u64> {
    pub original_commitments_root: [u8; 32],
    pub indirect_commitments: Vec<IndirectCommitment<N>>,
    pub direct_commitments: Vec<Commitment<N>>,
}

pub struct IndirectCommitment<const N: u64> {
    // Vec of diagonal commitments, each element contains 4 values over a 4-root in order to form a row
    pub diagonal_commitments: Vec<[Commitment<N>; 4]>,
    pub column_root: [u8; 32],
    // Vec of column commitments, each element contains 4 values over a 4-root, only the first one should be used in order to match the row from the diagonal commitment above
    pub column_commitments: Vec<Commitment<N>>,
}

fn derive_indirect_steps_count(max_degree: u64) -> Result<u32, anyhow::Error> {
    if max_degree == 0 {
        return Ok(0);
    }
    let mut degree = max_degree;
    let mut i = 0;
    while degree > DEGREE_THRESHOLD {
        if !degree.is_multiple_of(4) {
            return Err(anyhow!(
                "max_degree must be a divisible by 4 until degree threshold {DEGREE_THRESHOLD} is reached"
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
        result = 4 * result;
    }
    result
}
