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
        seed = Sha3_256::digest(seed).into();

        for chunk in seed.chunks(4) {
            let index = chunk_to_usize(chunk, modulus);
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
    chunk_to_usize(&digest[0..4], modulus)
}

fn chunk_to_usize(chunk: &[u8], modulus: usize) -> usize {
    let mut le_bytes = [0u8; 4];
    le_bytes.copy_from_slice(chunk);
    u32::from_le_bytes(le_bytes) as usize % modulus
}

#[cfg(test)]
mod prf_tests {
    use std::hash::{BuildHasherDefault, DefaultHasher};

    use super::*;

    #[test]
    fn test_prf_many_indices() {
        let seed: [u8; 32] = rand::random();

        let modulus = 100;
        let excluded_indices = HashSet::from_iter(0..5);

        let target_length = rand::random_range(15..25);

        let indices =
            pseudo_random_select_units_indices(&seed, target_length, modulus, &excluded_indices)
                .unwrap();
        assert!(indices.iter().any(|i| !excluded_indices.contains(i)));
        assert!(indices.len() == target_length);

        let selected_indices: HashSet<_, BuildHasherDefault<DefaultHasher>> =
            HashSet::from_iter(indices.iter());
        assert!(selected_indices.len() == target_length);
    }
}
