use std::collections::HashSet;

use anyhow::anyhow;
use sha3::{Digest, Sha3_256};

const MAX_ITERATIONS: usize = 50_000;

pub fn pseudo_random_select_units_indices(
    seed: &[u8; 32],
    target_length: usize,
    modulus: usize,
    excluded_indices: &HashSet<usize>,
) -> Result<Vec<usize>, anyhow::Error> {
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

            if !excluded_indices.contains(&index) && !selected_indices.contains(&index) {
                selected_indices.insert(index);
                indices.push(index);
            }
        }
        counter += 1;

        if counter > MAX_ITERATIONS {
            return Err(anyhow!(
                "exceeded maximum iterations for pseudo-random selection"
            ));
        }
    }
    indices.resize(target_length, 0);
    Ok(indices)
}

pub fn pseudo_random_select_unit_index(seed: &[u8; 32], modulus: usize) -> usize {
    let digest: [u8; 32] = Sha3_256::digest(seed).into();
    let mut le_bytes = [0u8; 4];
    le_bytes.clone_from_slice(&digest[0..4]);
    u32::from_le_bytes(le_bytes) as usize % modulus
}
