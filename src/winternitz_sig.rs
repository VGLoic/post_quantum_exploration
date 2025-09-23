#[cfg(test)]
mod test {
    use anyhow::anyhow;
    use bitvec::{prelude::*, view::BitView};
    use fake::{Fake, Faker};
    use sha3::{Digest, Sha3_256};

    /*
     * # About Winternitz One Time Signatures (WOTS)
     *
     * Parameters:
     *  - `w`: integer of the form `2^x`,
     *  - secret key size: 32 bytes,
     *  - hashing algorithm: SHA3 256
     *
     * Signature process of `msg`:
     * 1. we decompose `msg` in bits,
     * 2. we split the bits by chunks of size `x`,
     * 3. we initialize a checksum value to 0,
     * 4. for each chunk:
     *  4.a. we remap the bits into an integer, this integer is lower than `w = 2^x` by design, we call this mapped value `element`,
     *  4.b. we randomly generate a secret key,
     *  4.c. we generate the signature by hashing `element` times the secret key,
     *  4.d. we generate the `public key` associated to the secret key by hashing the secret key `w` times,
     *  4.e. we increment the checksum value with `w - 1 - element`,
     * 5. we repeat the steps 2. to 4. for the checksum, hence obtaining additional signatures with their secret keys and public keys,
     * 6. the output `sig` is:
     *  - the list of the signatures of each elements,
     *  - the list of the signatures of the checksum.
     *
     *
     * Verification process of `msg` and `sig`:
     * 1. we decompose `msg` in bits,
     * 2. we split the bits by chunks of size `x`,
     * 3. we initialize a checksum value to 0,
     * 4. for each chunk:
     *  4.a. we remap the bits into an integer, this integer is lower than `w = 2^x` by design, we call this mapped value `element`,
     *  4.b. we hash the associated input element signature `w - element` times,
     *  4.c. element signature is valid if obtained hash is equal to associated input element public key,
     *  4.d. we increment the checksum value with `w - 1 - element`,
     * 5. we repeat the steps 2. to 4. for the checksum in order to verify the checksums signatures.
     *
     *
     * Points to consider:
     * 1. In terms of size, given `w = 2^x`, we have `y = msg.len() / x` elements in a message, therefore:
     *  - we consume `y + checksum keys` secret keys,
     *  - the signature output is `y` element signatures (size of a hash), `y` public keys (size of a hash) and additional signatures and public keys for the checksum.
     * Examples with x = 4, w = 2^4 = 16, and msg of length 64,
     *  - `y = 64 / 4 = 16`
     *  - we got 16 element signatures,
     *  - checksum is at most `16 * w = 2^3 * 2^4 = 2^7` so we can decompose it in two groups of 4 bits. Therefore we have 2 additional elements to deal with,
     *  - in total: 18 (16 + 2) private keys consumed and as output: 18 public keys and 18 signatures
     * 2. Once we sign an element, e.g. `element = 5`, an attacker can re-use the keys (without knowing them) this to sign increment values of it by simply hashing the signature.
     *    Because of this, the `checksum` part was added: if an attacker increases an element, the checksum decreases and the attacker can not re-use the keys for the checksum to build a new valid one.
     *    The attacker would then need to generate new keys for the checksum part, but they would not be associated to the original signer.
     * 3. Because of the possible reusability of the keys, this signature scheme is meant to work by having unique keys per messages. There is a whole topic of key management for this:
     *  - need to generate a whole lot of keys,
     *  - need to keep track of which keys have been used or not.
     * It seems the more promising solutions are using Merkle trees to store the keys, hence having a single public key (the root hash) as the way to verify that the public key belongs to the signer.
     * => By doing so, we then need to replace the previous public keys in the output by the Merkle proofs of each key we use along the signature.
     *
     * Including Merkle tree leads to various possibilities with their tradeoff:
     * - generate one Merkle tree for each signature,
     * - generate one Merkle tree for multiple signatures. We can not really handle gigantic Merkle trees in memory like that (practical limit is depth of 16).
     *   So we introduce Multiple Merkle trees: a Merkle tree of Merkle trees where leaves are derived according to a pseudo random function and a seed. We can then generate the needed Merkle tree when needed but we don't need to store the whole structure.
     *   Navigating in these gigantic Merkle trees is a challenge so we got two variants:
     *      - Stateful hash-based signatures: we need to keep track of a counter in order to know which keys have been used or not. This is generally very complex to ensure so if people need this scheme, they even implement the key management at the hardware level,
     *      - Stateless hash-based signatures: the object to sign is used to derive which keys we will use.
     *        Example with `msg` of size `64`:
     *            - we want a tree that can sign all these messages: so we need 2^64 elements,
     *            - we use a multi Merkle tree of depth 4 (2^4 = 16 "leaves"), with Merkle tree of depth 16, hence we got a actual depth of 4 * 16 = 64, (and 2^64 leaves),
     *            - the `msg` is split in 4 chunks of 16 bits, each chunk is used to choose at each level which tree we use.
     *              E.g. first 16 bits for which tree at level 1, second 16 bits for which tree at level 2, third 16 bits for which tree at level 3, last 16 bits for which secret key at level 4.
     *        Related to the example, 64 bits signing is considered insecure, the minimum is 128 (still weak) and OK is 256. It leads to respectively MMS of depth 8 * 16 and 16 * 16 respectively. Therefore the signatures become bigger and bigger.
     *
     *
     * Related to this idea, the SPHINCS and SPHINCS+ are stateless signatures using this approach (variations as I understand).
     * A signature is between 8k and 49k bytes.
     */

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

    struct ElementSignature {
        sig: [u8; 32],
        public_key: [u8; 32],
    }

    struct Signature {
        element_signatures: Vec<ElementSignature>,
        checksum_signatures: Vec<ElementSignature>,
    }

    fn sign_element(
        element: u8,
        params: &WinternitzParameters,
    ) -> Result<ElementSignature, anyhow::Error> {
        if element >= params.w {
            return Err(anyhow!("message must be less than {}", params.w));
        }

        let secret_key: [u8; 32] = rand::random();

        // We hash `element` times in order to get the signature
        let mut sig = secret_key;
        for _ in 0..element {
            sig = Sha3_256::digest(sig).into();
        }

        // We then continue the hashing in order to get the public key defined as `secret key` hashes `w` times
        let mut public_key = sig;
        for _ in element..params.w {
            public_key = Sha3_256::digest(public_key).into();
        }

        Ok(ElementSignature { sig, public_key })
    }

    fn verify_element_signature(
        element: u8,
        signature: &ElementSignature,
        params: &WinternitzParameters,
    ) -> Result<(), anyhow::Error> {
        if element >= params.w {
            return Err(anyhow!("message must be less than {}", params.w));
        }

        let mut recovered_public_key = signature.sig;
        for _ in 0..(params.w - element) {
            recovered_public_key = Sha3_256::digest(recovered_public_key).into();
        }

        if recovered_public_key != signature.public_key {
            return Err(anyhow!("invalid signature"));
        }

        Ok(())
    }

    fn sign<T>(msg: T, params: &WinternitzParameters) -> Result<Signature, anyhow::Error>
    where
        T: BitView,
    {
        let mut checksum: u8 = 0;

        let mut element_signatures = vec![];
        // 64 / log2_w iterations -> 64 / log2_w signatures and public keys -> size
        // Example w = 16, log2_w = 4, hence we have 16 signatures and public keys of 256 bits each -> 16 * 256 = 4096 bits of signatures and public keys
        for element in to_padded_bits(msg, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
        {
            checksum += params.w - 1 - element;
            let element_signature = sign_element(element, params)?;
            element_signatures.push(element_signature);
        }

        let mut checksum_signatures = vec![];
        for element in to_padded_bits(checksum, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
        {
            let checksum_signature = sign_element(element, params)?;
            checksum_signatures.push(checksum_signature);
        }

        Ok(Signature {
            element_signatures,
            checksum_signatures,
        })
    }

    fn verify_signature<T>(
        msg: T,
        signature: Signature,
        params: &WinternitzParameters,
    ) -> Result<(), anyhow::Error>
    where
        T: BitView,
    {
        let mut checksum = 0u8;
        for (i, element) in to_padded_bits(msg, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
            .enumerate()
        {
            checksum += params.w - 1 - element;
            verify_element_signature(element, &signature.element_signatures[i], params)?;
        }

        for (i, element) in to_padded_bits(checksum, params)
            .chunks(params.log2_w.into())
            .map(map_bit_group_to_element)
            .enumerate()
        {
            verify_element_signature(element, &signature.checksum_signatures[i], params)?;
        }

        Ok(())
    }

    fn to_padded_bits<T>(
        msg: T,
        params: &WinternitzParameters,
    ) -> BitVec<<<T as BitView>::Store as BitStore>::Unalias>
    where
        T: BitView,
    {
        let mut bits = msg.view_bits::<Lsb0>().to_bitvec();
        let r = bits.len() / (params.log2_w as usize);
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
    fn test_signature_element() {
        let params = WinternitzParameters::new(4).unwrap();
        for msg in 0..params.w {
            let signature = sign_element(msg, &params).unwrap();
            assert!(verify_element_signature(msg, &signature, &params).is_ok());
        }
    }

    #[test]
    fn test_signature_larger_message() {
        let params: WinternitzParameters = WinternitzParameters::new(4).unwrap();
        let msg = Faker.fake::<u64>();
        assert!(verify_signature(msg, sign(msg, &params).unwrap(), &params).is_ok());
    }
}
