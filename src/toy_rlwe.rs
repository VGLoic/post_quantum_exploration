#[cfg(test)]
mod toy_rlwe_tests {

    /*
     * This is a toy example for Ring Learning With Errors (RLWE), not optimized for performance or size.
     * It illustrates the basic principles of RLWE-based key exchange using simple polynomials operations.
     *
     * The LWE (see toy_lwe.rs) protocol requires large vectors and matrices to be safe.
     * The RLWE comes as an attempt to decrease the size by replacing the vectors and matrices with polynomials on a ring, i.e. operations are done modulo a polynomial, e.g. x^n + 1.
     *
     * The protocol is a bit tricky as I read that the formal proof of difficulty no longer hold because we can't show equivalence between RLWE and the Lattice problems.
     */

    use crate::primefield::PrimeFieldElement;

    type FieldElement = PrimeFieldElement<1_000_000_007>;

    const RING_DEGREE: usize = 7;

    #[derive(Debug, Clone, PartialEq)]
    struct Polynomial {
        // coefficients are in ascending degrees, i.e. a0 + a_1 * x + a_2 * x^2 + ... + a_6 * x^6
        coefficients: [FieldElement; RING_DEGREE],
    }

    impl Polynomial {
        fn new(coefficients: [FieldElement; RING_DEGREE]) -> Self {
            Self { coefficients }
        }

        fn add(&self, other: &Self) -> Self {
            let mut result_coefficients = [FieldElement::new(0); RING_DEGREE];

            for (i, pair) in self
                .coefficients
                .iter()
                .zip(&other.coefficients)
                .enumerate()
            {
                result_coefficients[i] = pair.0.add(pair.1);
            }

            Self::new(result_coefficients)
        }

        fn mul(&self, other: &Self) -> Self {
            // We first compute the resulting coefficients, they can be at most of degree 2 * (RING_DEGREE - 1), hence 2 * RING_DEGREE - 1 coefficients
            let mut raw_coefficients = [FieldElement::new(0); 2 * RING_DEGREE - 1];
            for (i, a_i) in self.coefficients.iter().enumerate() {
                for (j, b_j) in other.coefficients.iter().enumerate() {
                    raw_coefficients[i + j] = raw_coefficients[i + j].add(&a_i.mul(b_j));
                }
            }

            let (coefficients, to_be_filtered) = raw_coefficients.split_at_mut(RING_DEGREE);

            for (i, coeff) in to_be_filtered.iter().enumerate() {
                // We subtract `coeff * x^i * (x^RING_DEGREE + 1)` to the full polynomial
                // The `RING_DEGREE + i` coefficient will then vanish
                // We only register the subtraction of the `i`th coefficient by `coeff`
                coefficients[i] = coefficients[i].add(&coeff.neg());
            }

            Self::new(coefficients.try_into().unwrap_or_else(|_| {
                unreachable!("coefficients is necessarily of length RING_DEGREE")
            }))
        }

        fn truncate_below(&self, degree_threshold: u8) -> Self {
            let mut coefficients = [FieldElement::new(0); RING_DEGREE];
            if degree_threshold as usize >= RING_DEGREE {
                return Self::new(coefficients);
            }
            coefficients[degree_threshold as usize..RING_DEGREE]
                .copy_from_slice(&self.coefficients[degree_threshold as usize..RING_DEGREE]);
            Self::new(coefficients)
        }

        fn random(degree: u8) -> Self {
            let mut coefficients = [FieldElement::new(0); RING_DEGREE];
            for coeff in coefficients
                .iter_mut()
                .take((degree as usize + 1).min(RING_DEGREE))
            {
                *coeff = FieldElement::new(rand::random());
            }
            Self::new(coefficients)
        }
    }

    impl std::fmt::Display for Polynomial {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let terms: Vec<String> = self
                .coefficients
                .iter()
                .enumerate()
                .filter(|(_, coeff)| coeff.inner() != 0)
                .map(|(i, coeff)| {
                    if i == 0 {
                        format!("{}", coeff)
                    } else if i == 1 {
                        format!("{}x", coeff)
                    } else {
                        format!("{}x^{}", coeff, i)
                    }
                })
                .collect();
            write!(f, "{}", terms.join(" + "))
        }
    }

    #[test]
    fn test_rlwe_dh() {
        let a = Polynomial::new([
            0.into(),
            0.into(),
            0.into(),
            0.into(),
            FieldElement::new(rand::random()),
            FieldElement::new(rand::random()),
            FieldElement::new(rand::random()),
        ]);

        let s_alice = Polynomial::random(2);
        let e_alice = Polynomial::random(2);
        let p_alice = a.mul(&s_alice).add(&e_alice);

        let s_bob = Polynomial::random(2);
        let e_bob = Polynomial::random(2);
        let p_bob = a.mul(&s_bob).add(&e_bob);

        let k_alice = p_bob.mul(&s_alice).truncate_below(5);
        let k_bob = p_alice.mul(&s_bob).truncate_below(5);

        assert_eq!(k_alice, k_bob);
    }

    #[test]
    fn test_polynomial_addition() {
        let a = FieldElement::new(rand::random());
        let b = FieldElement::new(rand::random());
        let p = Polynomial::new([a, a, a, a, a, a, a]);
        let q = Polynomial::new([b, b, b, b, b, b, b]);
        let a_plus_b = a.add(&b);
        assert_eq!(
            p.add(&q),
            Polynomial::new([
                a_plus_b, a_plus_b, a_plus_b, a_plus_b, a_plus_b, a_plus_b, a_plus_b,
            ])
        );
    }

    #[test]
    fn test_polynomial_multiplication() {
        let a = FieldElement::new(5);
        let b = FieldElement::new(6);
        let zero = FieldElement::new(0);
        let p = Polynomial::new([a, zero, zero, zero, zero, zero, zero]);
        let q = Polynomial::new([b, zero, zero, zero, zero, zero, zero]);
        assert_eq!(
            p.mul(&q),
            Polynomial::new([a.mul(&b), zero, zero, zero, zero, zero, zero,])
        );

        // P(x) = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
        // Q(x) = 2 + 4x + 8x^2 + 16x^3 + 32x^4 + 64x^5 + 128x^6
        // (P * Q)(x) = 2 + 6x + 14x^2 + 30x^3 + 62x^4 + 126x^5 + 254x^6
        //           + 252x^7 + 248x^8 + 240x^9 + 224x^{10} + 192x^{11} + 128x^{12}
        // After reduction mod x^7 + 1:
        // (P * Q)(x) ≡ (2 - 252) + (6 - 248)x + (14 - 240)x^2 + (30 - 224)x^3
        //           + (62 - 192)x^4 + (126 - 128)x^5 + (254)x^6
        //         = -250 + (-242)x + (-226)x^2 + (-194)x^3 + (-130)x^4 + (-2)x^5 + 254x^6
        // (all coefficients mod p, if using a prime field)
        // For p = 1_000_000_007:
        // (P * Q)(x) ≡ 999999757 + 999999765x + 999999781x^2 + 999999813x^3 + 999999877x^4 + 1000000005x^5 + 254x^6
        let one = FieldElement::new(1);
        let p = Polynomial::new([one, one, one, one, one, one, one]);
        let q = Polynomial::new([
            2.into(),
            4.into(),
            8.into(),
            16.into(),
            32.into(),
            64.into(),
            128.into(),
        ]);
        assert_eq!(
            p.mul(&q),
            Polynomial::new([
                FieldElement::new(999_999_757),   // -250 mod 1_000_000_007
                FieldElement::new(999_999_765),   // -242 mod 1_000_000_007
                FieldElement::new(999_999_781),   // -226 mod 1_000_000_007
                FieldElement::new(999_999_813),   // -194 mod 1_000_000_007
                FieldElement::new(999_999_877),   // -130 mod 1_000_000_007
                FieldElement::new(1_000_000_005), // -2 mod 1_000_000_007
                FieldElement::new(254),
            ])
        );
    }
}
