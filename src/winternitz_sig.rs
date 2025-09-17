#[cfg(test)]
mod test {
    use anyhow::anyhow;
    use bitvec::{prelude::*, view::BitView};
    use fake::{Fake, Faker};
    use sha3::{Digest, Sha3_256};

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

        let mut secret_key: [u8; 32] = rand::random();
        // The first byte of the secret key is set to the message
        secret_key[0] = element;

        // We hash `element` times in order to get the signature
        let mut sig = secret_key;
        for _ in 0..element {
            sig = Sha3_256::digest(sig).into();
        }

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
            checksum += element;
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
            checksum += element;
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
