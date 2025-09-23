/*
 * Toy example of stateless Winternitz signatures using
 * - message size of 8 bits,
 * - Winternitz parameters: w = 16 = 2^4
 * - we use a Merkle tree of depth 8 to handle the 2^8 possible keys,
 *
 *
 * We end up with groups of 4 bits, hence y = 8 / 4 = 2 and checksum is at most 2 * 16 = 32 = 2^6
 * Therefore: 4 (2 + 2) signatures and public keys
 */

#[cfg(test)]
mod stateless_winternitz_tests {
    use anyhow::anyhow;
    use bitvec::{prelude::*, view::BitView};
    use sha3::{Digest, Sha3_256};

    use crate::merkletree_v2::{MerkleTreeV2, verify_proof};

    #[derive(Debug, Clone)]
    struct KeyPack {
        secret_keys: [[u8; 32]; 4],
        packed_public_keys: Vec<u8>,
    }

    impl KeyPack {
        fn new(hashing_round: u8) -> Self {
            let mut secret_keys: [[u8; 32]; 4] = [[0u8; 32]; 4];
            let mut packed_public_keys = vec![];
            for item in &mut secret_keys {
                let secret_key: [u8; 32] = rand::random();
                let mut public_key = Sha3_256::digest(secret_key);
                for _ in 1..hashing_round {
                    public_key = Sha3_256::digest(public_key);
                }
                for u in public_key.as_slice() {
                    packed_public_keys.push(*u);
                }
                *item = secret_key;
            }
            Self {
                secret_keys,
                packed_public_keys,
            }
        }
    }

    impl Default for KeyPack {
        fn default() -> Self {
            Self {
                secret_keys: [[0; 32]; 4],
                packed_public_keys: vec![0u8; 128],
            }
        }
    }

    impl AsRef<[u8]> for KeyPack {
        fn as_ref(&self) -> &[u8] {
            &self.packed_public_keys
        }
    }

    struct WinternitzParameters {
        w: u8,
        log2_w: u8,
    }

    impl WinternitzParameters {
        fn new(log2_w: u8) -> Result<Self, anyhow::Error> {
            let w = (2u8)
                .checked_pow(log2_w.into())
                .ok_or(anyhow!("invalid parameters, got a too big w"))?;
            Ok(WinternitzParameters { w, log2_w })
        }
    }

    struct Signature {
        element_signatures: Vec<[u8; 32]>,
        checksum_signatures: Vec<[u8; 32]>,
        merkle_proof: Vec<(bool, [u8; 32])>,
    }

    fn sign(
        msg: u8,
        params: &WinternitzParameters,
        tree: &MerkleTreeV2<KeyPack>,
    ) -> Result<Signature, anyhow::Error> {
        let proven_key_pack = tree.get(&msg_to_selector(msg))?;
        let mut checksum: u8 = 0;

        let mut secret_keys_iterator = proven_key_pack.value.secret_keys.iter();

        let mut element_signatures = vec![];
        for element in to_padded_bits(msg, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
        {
            checksum += params.w - 1 - element;
            let element_signature = sign_element(
                element,
                secret_keys_iterator
                    .next()
                    .ok_or(anyhow!("unexpected end of iterator"))?,
                params,
            )?;
            element_signatures.push(element_signature);
        }

        let mut checksum_signatures = vec![];
        for element in to_padded_bits(checksum, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
        {
            let checksum_signature = sign_element(
                element,
                secret_keys_iterator
                    .next()
                    .ok_or(anyhow!("unexpected end of iterator"))?,
                params,
            )?;
            checksum_signatures.push(checksum_signature);
        }

        Ok(Signature {
            element_signatures,
            checksum_signatures,
            merkle_proof: proven_key_pack
                .proof
                .into_iter()
                .map(|v| (v.0, v.1.to_owned()))
                .collect(),
        })
    }

    fn verify_signature(
        msg: u8,
        signature: Signature,
        root_hash: &[u8; 32],
        params: &WinternitzParameters,
    ) -> Result<(), anyhow::Error> {
        let mut checksum = 0u8;
        let mut packed_public_keys = vec![];
        for (i, element) in to_padded_bits(msg, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
            .enumerate()
        {
            checksum += params.w - 1 - element;

            let mut recovered_public_key = signature.element_signatures[i];
            for _ in 0..(params.w - element) {
                recovered_public_key = Sha3_256::digest(recovered_public_key).into();
            }
            for u in recovered_public_key {
                packed_public_keys.push(u);
            }
        }

        for (i, element) in to_padded_bits(checksum, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
            .enumerate()
        {
            let mut recovered_public_key = signature.checksum_signatures[i];
            for _ in 0..(params.w - element) {
                recovered_public_key = Sha3_256::digest(recovered_public_key).into();
            }
            for u in recovered_public_key {
                packed_public_keys.push(u);
            }
        }

        verify_proof(root_hash, &packed_public_keys, &signature.merkle_proof)
    }

    fn sign_element(
        element: u8,
        secret_key: &[u8; 32],
        params: &WinternitzParameters,
    ) -> Result<[u8; 32], anyhow::Error> {
        if element >= params.w {
            return Err(anyhow!("message must be less than {}", params.w));
        }

        // We hash `element` times in order to get the signature
        let mut sig = secret_key.to_owned();
        for _ in 0..element {
            sig = Sha3_256::digest(sig).into();
        }

        Ok(sig)
    }

    fn msg_to_selector(msg: u8) -> [bool; 8] {
        let vec_selector: Vec<bool> = msg.view_bits::<Lsb0>().to_bitvec().into_iter().collect();
        let mut selector = [false; 8];
        selector.clone_from_slice(&vec_selector);
        selector
    }

    fn to_padded_bits(
        msg: u8,
        params: &WinternitzParameters,
    ) -> BitVec<<<u8 as BitView>::Store as BitStore>::Unalias> {
        let mut bits = msg.view_bits::<Lsb0>().to_bitvec();

        let r = bits.len() % (params.log2_w as usize);
        if r != 0 {
            bits.resize(bits.len() + (params.log2_w as usize) - r, false);
        }

        bits
    }

    fn map_bit_group_to_element<T>(bit_group: &BitSlice<T, Lsb0>) -> u8
    where
        T: bitvec::store::BitStore,
    {
        let mut acc: u8 = 0;
        for (i, bit) in bit_group.iter().enumerate() {
            if bit == true {
                acc += (2u8).pow(i.try_into().unwrap())
            }
        }
        acc
    }

    #[test]
    fn test_stateless_winternitz_sig() {
        let params = WinternitzParameters::new(4).unwrap();

        let mut values = vec![];
        for _ in 0..(2usize.pow(8)) {
            values.push(KeyPack::new(params.w));
        }
        let tree = MerkleTreeV2::new(8, &values).unwrap();

        let msg: u8 = rand::random();
        let signature = sign(msg, &params, &tree).unwrap();
        verify_signature(msg, signature, tree.root_hash(), &params).unwrap();
    }
}
