use std::collections::HashSet;

use sha3::{Digest, Sha3_256};

pub fn pseudo_random_select_units_indices(
    seed: &[u8; 32],
    target_length: usize,
    modulus: usize,
) -> Vec<usize> {
    let mut indices = vec![];
    let mut selected_indices: HashSet<usize> = HashSet::new();
    let mut counter = 0;
    let mut seed: [u8; 32] = *seed;
    while indices.len() < target_length {
        seed = Sha3_256::digest(seed).into(); // 32 bytes
        for chunk in seed.chunks(4) {
            let mut le_bytes = [0u8; 4];
            le_bytes.clone_from_slice(chunk);
            let index = u32::from_le_bytes(le_bytes) as usize % modulus;

            if !selected_indices.contains(&index) {
                selected_indices.insert(index);
                indices.push(index);
            }
        }
        counter += 1;

        if counter > 10_000 {
            panic!("oula");
        }
    }
    indices.resize(target_length, 0);
    indices
}

pub fn pseudo_random_select_indirect_proof_indices(
    seed: &[u8; 32],
    row_count: usize,
    modulus: usize,
) -> (usize, Vec<usize>) {
    let mut indices = pseudo_random_select_units_indices(seed, row_count + 1, modulus);
    let column_index = indices.pop().unwrap();
    (column_index, indices)
}
