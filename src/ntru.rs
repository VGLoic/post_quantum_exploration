#[cfg(test)]
mod ntru_tests {

    use crate::centered_modular_arithmetic::*;

    /*
     * This is a simple implementation of the NTRU encryption scheme
     * See https://en.wikipedia.org/wiki/NTRUEncrypt for more details
     *
     * This implementation is not optimized for performance or security, and is only meant for educational purposes.
     * In particular, the random polynomial generation is not secure, and the parameters are not chosen for security.
     *
     * This scheme is based on polynomial arithmetic in the ring R = Z[x] / (x^N + 1) where N is prime.
     *
     * We first generate polynomials with low amplitude coefficients in {-1, 0, 1} for the private keys or the messages to encrypt.
     * The public keys or the cyphertexts are then polynomials with larger coefficients in a larger modulus.
     * Knowledge of the private key allows to reduce the coefficients back to the low amplitude range and recover the message.
     *
     * In more details:
     * # Key generation
     * - Generate two polynomials `f` and `g` with coefficients in {-1, 0, 1},
     * - Derive `f_p` and `f_q`, the inverses of `f` modulo `p` and `q` respectively,
     * - The public key is `h = p * g * f_q (mod q)`,
     * - The private key is `(f, f_p)`.
     *
     * # Encryption
     * - Generate a random polynomial `r` with coefficients in {-1, 0, 1},
     * - The cyphertext is `c = public_key * r + m (mod q)`, where `m` is the message polynomial with coefficients in {-1, 0, 1}.
     *
     * # Decryption
     * - Compute `v = c * f (mod q) = p * r * g + f * m (mod q)`,
     * - Compute `v * f_p (mod p) = m (mod p)`, and recover `m`.
     *
     * # Note 1
     * Because of the ring structure, the multiplication is done modu `x^N + 1` and the degree of the polynomials is bounded by `N - 1`.
     *
     * # Note 2
     * The operations are done modulo `p` or `q`, and the coefficients are always reduced to the centered range `[-(n/2), (n/2)]`.
     * This centered reduction is crucial to ensure the correctness of the decryption. Otherwise, the coefficients could wrap around and lead to incorrect results.
     * In essence, the cyphertext computation is done modulo `q`, but the coefficients are expected to be in the range `[-(q/2), (q/2)]` to ensure that when multiplied by `f` and reduced modulo `p`, the original message `m` can be correctly recovered.
     */

    // NTRU parameters
    //  - N: polynomial degree bound
    //  - Q: large modulus,
    //  - P: small modulus.
    // N is prime, P and Q are coprime
    const N: usize = 7;
    const Q: u32 = 2048;
    const P: u32 = 3;

    #[derive(Debug)]
    struct NtruKeyPair {
        pub public_key: Polynomial,
        pub private_f: Polynomial,
        pub private_f_p: Polynomial,
    }

    impl NtruKeyPair {
        /// Generate a new NTRU key pair
        /// - The public key is `h = p * g * f_q (mod q)`, where `f_q` is the inverse of `f mod q`
        /// - The private key is `(f, f_p)`, where `f_p` is the inverse of `f mod p`
        ///
        /// The polynomials `f` and `g` are randomly generated with coefficients in {-1, 0, 1}
        ///
        /// We ensure that `f` is invertible mod `p` and mod `q`, and that `g` is not the zero polynomial.
        fn generate() -> Self {
            let mut g = Polynomial::generate(P);
            while g.is_zero() {
                g = Polynomial::generate(P);
            }

            let mut f = Polynomial::default();
            let mut f_p = Polynomial::default();
            let mut f_q = Polynomial::default();
            while f_p.is_zero() && f_q.is_zero() {
                f = Polynomial::generate(P);
                if let Some(p) = f.inv(P)
                    && let Some(q) = f.inv(Q)
                {
                    f_p = p;
                    f_q = q;
                }
            }

            let public_key = g.mul(&f_q, Q).scalar_mul(P.into(), Q);

            Self {
                public_key,
                private_f: f,
                private_f_p: f_p,
            }
        }

        /// Decrypt a cyphertext using a private key and outputs the decrypted message
        ///
        /// Let us consider the definition of the cyphertext
        /// ```
        /// c = public_key * r + m (mod q) = p * r * g * f_q + m (mod q)
        /// ```
        /// We multiply by `c` by `f`:
        /// ```
        /// v = c * f (mod q) = p * r * g * f * f_q + f * m (mod q) = p * r * g + f * m (mod q)
        /// ```
        /// We then multiply by f_p mod p
        /// ```
        /// v * f_p (mod p) = p * r * g * f_p + f * f_p * m (mod p) = m (mod p)
        /// ```
        /// And `m` is recovered.
        ///
        ///
        /// # Arguments
        /// * `cyphertext` - the cyphertext to decrypt
        fn decrypt(&self, cyphertext: &Polynomial) -> Polynomial {
            let v = cyphertext.mul(&self.private_f, Q);
            v.mul(&self.private_f_p, P)
        }
    }

    /// Encrypt a message using the recipient public key. The cyphertext is returned.
    ///
    /// The cyphertext is defined as:
    /// ```
    /// c = public_key * r + m (mod q)
    /// ```
    /// Where `r` is a randomly generated polynomial.
    ///
    /// # Arguments
    /// * `public_key` - the recipient public key,
    /// * `message` - the message to encrypt, encoded as a polynomial
    fn encrypt(public_key: &Polynomial, message: &Polynomial) -> Polynomial {
        let r = Polynomial::generate(P);

        r.mul(public_key, Q).add(message, Q)
    }

    #[derive(Debug, Clone, PartialEq, Default)]
    struct Polynomial {
        coefficients: [i64; N],
    }

    impl std::fmt::Display for Polynomial {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let terms: Vec<String> = self
                .coefficients
                .iter()
                .enumerate()
                .map(|(i, &coeff)| {
                    if coeff == 0 {
                        "".to_string()
                    } else if i == 0 {
                        format!("{}", coeff)
                    } else if i == 1 {
                        format!("{}x", coeff)
                    } else {
                        format!("{}x^{}", coeff, i)
                    }
                })
                .filter(|d| !d.is_empty())
                .collect();
            if terms.is_empty() {
                return write!(f, "0");
            }
            write!(f, "{}", terms.join(" + "))
        }
    }

    impl Polynomial {
        /// Create a new polynomial from the given coefficients
        /// # Arguments
        /// * `coefficients` - the coefficients of the polynomial, in increasing order of degree
        /// # Example
        /// ```
        /// let p = Polynomial::new([1, 2, 3, 4, 5, 6, 7]); // 1 + 2x + 3x^2 + 4x^3 + 5x^4 + 6x^5 + 7x^6
        /// ```
        fn new(coefficients: [i64; N]) -> Self {
            Self { coefficients }
        }

        /// Generate a random polynomial with coefficients in the range `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `n` - the modulus to use for the coefficients
        /// # Example
        /// ```
        /// let p = Polynomial::generate(13); // Random polynomial with coefficients in [-6, 6]
        /// ```
        fn generate(n: u32) -> Self {
            let mut coefficients = [0i64; N];
            for coeff in coefficients.iter_mut() {
                *coeff = to_centered_coordinates(rand::random(), n);
            }
            Self::new(coefficients)
        }

        /// Add two polynomials modulo `n`, resulting polynomial has coefficients in `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `other` - the other polynomial to add
        /// * `n` - the modulus to use for the addition
        fn add(&self, other: &Self, n: u32) -> Self {
            let mut result_coefficients = [0i64; N];

            for (i, pair) in self
                .coefficients
                .iter()
                .zip(&other.coefficients)
                .enumerate()
            {
                result_coefficients[i] = modulo_add(*pair.0, *pair.1, n);
            }

            Self::new(result_coefficients)
        }

        /// Negate a polynomial modulo `n`, resulting polynomial has coefficients in `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `n` - the modulus to use for the negation
        fn neg(&self, n: u32) -> Self {
            let mut coefficients = [0i64; N];

            for (i, coeff) in self.coefficients.iter().enumerate() {
                coefficients[i] = modulo_neg(*coeff, n);
            }

            Self::new(coefficients)
        }

        /// Check if the polynomial is the zero polynomial
        fn is_zero(&self) -> bool {
            !self.coefficients.iter().any(|v| *v != 0)
        }

        /// Return the degree of the polynomial
        fn degree(&self) -> usize {
            for i in (0..N).rev() {
                if self.coefficients[i] != 0 {
                    return i;
                }
            }
            0
        }

        /// Multiply two polynomials modulo `n` and modulo `x^N + 1`, resulting polynomial has coefficients in `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `other` - the other polynomial to multiply
        /// * `n` - the modulus to use for the multiplication
        fn mul(&self, other: &Self, n: u32) -> Self {
            let mut raw_coefficients = [0i64; 2 * N - 1];
            for (i, a_i) in self.coefficients.iter().enumerate() {
                for (j, b_j) in other.coefficients.iter().enumerate() {
                    raw_coefficients[i + j] =
                        modulo_add(raw_coefficients[i + j], modulo_mul(*a_i, *b_j, n), n);
                }
            }

            let (coefficients, to_be_filtered) = raw_coefficients.split_at_mut(N);

            for (i, coeff) in to_be_filtered.iter().enumerate() {
                // We subtract `coeff * x^i * (x^N + 1)` to the full polynomial
                // The `N + i` coefficient will then vanish
                // We only register the subtraction of the `i`th coefficient by `coeff`
                coefficients[i] = modulo_add(coefficients[i], modulo_neg(*coeff, n), n);
            }

            Self::new(
                coefficients
                    .try_into()
                    .unwrap_or_else(|_| unreachable!("coefficients is necessarily of length N")),
            )
        }

        /// Multiply a polynomial by a scalar modulo `n`, resulting polynomial has coefficients in `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `scalar` - the scalar to multiply by
        /// * `n` - the modulus to use for the multiplication
        fn scalar_mul(&self, scalar: i64, n: u32) -> Self {
            let mut coefficients = [0i64; N];
            for (i, &coeff) in self.coefficients.iter().enumerate() {
                coefficients[i] = modulo_mul(coeff, scalar, n);
            }
            Polynomial::new(coefficients)
        }

        /// Divide two polynomials modulo `n`, resulting (quotient, remainder) tuple of polynomials has coefficients in `[-(n/2), (n/2)]`
        /// # Arguments
        /// * `other` - the other polynomial to divide by
        /// * `n` - the modulus to use for the division
        fn div(&self, other: &Self, n: u32) -> Option<(Self, Self)> {
            let other_degree = other.degree();
            let self_degree = self.degree();
            if self.is_zero() || other.is_zero() || other_degree > self_degree {
                return Some((Polynomial::new([0i64; N]), self.clone()));
            }

            // Starting from the highest degree term of self, we will subtract multiples of other until we reach the degree of other
            // e.g. if self has degree 6 and other has degree 3, we will do the iterationss for degrees 6, 5, 4 and 3.

            let mut remainder = self.clone();
            let mut quotient_coefficients = [0i64; N];

            let other_leading_coefficient = other.coefficients[other_degree];

            for i in (other_degree..=self_degree).rev() {
                let other_leading_coefficient_inv = modulo_inv(other_leading_coefficient, n)?;
                let quotient_coefficient =
                    modulo_mul(remainder.coefficients[i], other_leading_coefficient_inv, n);

                quotient_coefficients[i - other_degree] = quotient_coefficient;

                let mut to_be_subtracted_coeff = [0i64; N];
                to_be_subtracted_coeff[i - other_degree] = quotient_coefficient;
                let to_be_subtracted = Polynomial::new(to_be_subtracted_coeff).mul(other, n).neg(n);

                remainder = remainder.add(&to_be_subtracted, n);
            }

            Some((Polynomial::new(quotient_coefficients), remainder))
        }

        /// Compute the inverse of the polynomial modulo `n` and modulo `x^N + 1`, resulting polynomial has coefficients in `[-(n/2), (n/2)]`
        /// If the polynomial is not invertible, returns None
        /// # Arguments
        /// * `n` - the modulus to use for the inversion
        fn inv(&self, n: u32) -> Option<Polynomial> {
            if self.is_zero() {
                return None;
            }
            let self_degree = self.degree();

            if self_degree == 0 {
                let inv_coeff = modulo_inv(self.coefficients[0], n)?;
                return Some(Polynomial::new([inv_coeff, 0, 0, 0, 0, 0, 0]));
            }

            // We can't handle `R(x) = x^N + 1` in our struct so we will perform a first division step of `R` by `self` in order to remove the x^N coefficient
            // We end up with larger = R - correction * self
            let (corrected, correction) = {
                let self_leading_coefficient = self.coefficients[self_degree];
                let quotient_coefficient = modulo_inv(self_leading_coefficient, n)?;

                let mut to_be_subtracted_coeff = [0i64; N];
                to_be_subtracted_coeff[N - self_degree] = quotient_coefficient;
                let correction = Polynomial::new(to_be_subtracted_coeff);

                // Every coefficient is multiplied by `quotient_coefficient`, increase its order of `N - self_degree` and is negated, the highest one is ignored
                let degree_increase = N - self_degree;
                let mut corrected_coefficients = [0i64; N];
                for (i, coeff) in self.coefficients.iter().enumerate() {
                    if i + degree_increase < N {
                        corrected_coefficients[i + degree_increase] =
                            modulo_neg(modulo_mul(*coeff, quotient_coefficient, n), n);
                    }
                }
                // Coming from the ring
                corrected_coefficients[0] += 1;

                (Polynomial::new(corrected_coefficients), correction)
            };

            let (larger, mut smaller, has_been_swapped) = if corrected.degree() > self_degree {
                (&corrected, self.clone(), false)
            } else {
                (self, corrected, true)
            };

            let (mut q, mut r) = larger.div(&smaller, n)?;

            if r.is_zero() {
                if smaller.degree() != 0 {
                    return None;
                }
                // We have R - correction * self = b
                // Hence the inversion is `- (1/b) * correction`
                let factor = modulo_neg(modulo_inv(smaller.coefficients[0], n)?, n);
                return Some(correction.mul(&Polynomial::new([factor, 0, 0, 0, 0, 0, 0]), n));
            }

            let mut quotients = vec![];

            // While the remainder is not zero, we register the quotient and we re-iterate with `(smaller, r = larger % smaller)`
            while !r.is_zero() {
                quotients.push(q);
                let larger = smaller;
                smaller = r;
                (q, r) = larger.div(&smaller, n)?;
            }

            // There is no inverse if the GCD of self with the ring polynomial is not one (or a constant)
            if smaller.degree() != 0 {
                return None;
            }

            // Now we can rebuild the inverse by reversing the loop above
            // At this stage we have larger_k = smaller_k * q_k + 0
            // Since smaller_k = r_(k-1) and smaller_k is of degree 0, we have at the last iteration):
            // The last recorded iteration k gives us:
            // larger_k = smaller_k * q_k * r_k
            // Since `smaller_k = r_(k-1)` and last smaller is of degree 0, we have `r_k` of degree 0, therefore we obtain the equivalence with:
            // 1 = (1/r_k) * (larger_k - smaller_k * q_k)
            // Let us now define `alpha_k` and `beta_k` as 1 = alpha_k * larger_k + beta_k * smaller_k
            // Injecting the relations `smaller_k = r_(k-1) = larger_(k-1) - q_(k-1) * smaller_(k-1)` and `larger_k = smaller_(k-1)`,
            // We obtain the suite definitions for `alpha_k` and `beta_k`:
            // ```
            // alpha_(i-1) = beta_i; alpha_k = (1/r_k),
            // beta_(i-1) = alpha_i - q_(i-1) * beta_i; beta_k = - (1/r_k) * q_k
            // ```

            let last_quotient = quotients
                .pop()
                .unwrap_or_else(|| unreachable!("quotient is necessarily filled with one value"));

            let r_k_inv = modulo_inv(smaller.coefficients[0], n)?;

            let mut alpha = Polynomial::new([r_k_inv, 0, 0, 0, 0, 0, 0]);
            let mut beta = last_quotient.mul(&alpha.neg(n), n);

            for quotient in quotients.iter().rev() {
                let new_beta = alpha.add(&quotient.mul(&beta, n).neg(n), n);
                alpha = beta;
                beta = new_beta;
            }

            // We adjust with respect to the previous correction:
            //  - If there has been a swap at the start, then larger = self and smaller = R - correction * self
            //    So we have 1 = alpha * self + beta * (R - correction * self) <=> 1 = beta * R + (alpha - beta * correction) * self
            //  - Otherwise larger = R - correction * self and smaller = self
            //    So we have 1 = alpha * (R - correction * self) + beta * self <=> 1 = alpha * R + (beta - alpha * correction) * self
            if has_been_swapped {
                Some(alpha.add(&beta.mul(&correction, n).neg(n), n))
            } else {
                Some(beta.add(&alpha.mul(&correction, n).neg(n), n))
            }
        }
    }

    #[test]
    fn test_ntru_encryption() {
        let message = Polynomial::generate(P);
        let key_pair = NtruKeyPair::generate();

        let cyphertext = encrypt(&key_pair.public_key, &message);
        assert_eq!(key_pair.decrypt(&cyphertext), message);
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new([1, 2, 3, 4, 5, 6, 7]);
        let p2 = Polynomial::new([7, 6, 5, 4, 3, 2, 1]);
        let expected = Polynomial::new([-5, -5, -5, -5, -5, -5, -5]);
        assert_eq!(p1.add(&p2, 13), expected);

        let p3 = Polynomial::new([12, 12, 12, 12, 12, 12, 12]);
        let p4 = Polynomial::new([1, 1, 1, 1, 1, 1, 1]);
        let expected2 = Polynomial::new([0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(p3.add(&p4, 13), expected2);

        let p3 = Polynomial::new([12, 12, 12, 12, 12, 12, 12]);
        let p4 = Polynomial::new([1, 1, 1, 1, 1, 1, 1]);
        let expected2 = Polynomial::new([0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(p3.add(&p4, 13), expected2);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new([1, 1, 1, 1, 1, 1, 1]);
        let p2 = Polynomial::new([2, 4, 8, 16, 32, 64, 128]);
        let expected = Polynomial::new([-250, -242, -226, -194, -130, -2, 254]);
        assert_eq!(p1.mul(&p2, 2048), expected);
    }

    #[test]
    fn test_polynomial_division() {
        let p1 = [0, 0, 0, 1, 0, 0, 1]; // x^6 + x^3
        let p2 = [1, 0, 0, 1, 0, 0, 0]; // x^3 + 1
        let (quotient, remainder) = Polynomial::new(p1).div(&Polynomial::new(p2), 13).unwrap();
        assert_eq!(quotient, Polynomial::new([0, 0, 0, 1, 0, 0, 0])); // x^3
        assert_eq!(remainder, Polynomial::new([0, 0, 0, 0, 0, 0, 0])); // 0

        let p1 = [1, 2, 0, 0, 0, 0, 1]; // x^6 + 2x + 1
        let p2 = [1, 0, 0, 1, 0, 0, 0]; // x^3 + 1
        // x^6 + 2x + 1 = x^3 * (x^3 + 1) -x^3 + 2x + 1 = (x^3 - 1) * (x^3 + 1) + 2x + 2
        let (quotient, remainder) = Polynomial::new(p1).div(&Polynomial::new(p2), 13).unwrap();
        assert_eq!(quotient, Polynomial::new([-1, 0, 0, 1, 0, 0, 0])); // x^3 - 1
        assert_eq!(remainder, Polynomial::new([2, 2, 0, 0, 0, 0, 0])); // 2x + 2
    }

    #[test]
    fn test_polynomial_inv() {
        let one = Polynomial::new([1, 0, 0, 0, 0, 0, 0]);
        for p in [
            one.clone(),
            Polynomial::new([8, 0, 0, 0, 0, 0, 0]), // 1
            Polynomial::new([2, 0, 0, 0, 0, 0, 0]), // 2
            Polynomial::new([0, 1, 0, 0, 0, 0, 0]), // x
            Polynomial::new([0, 8, 0, 0, 0, 0, 0]), // 8x
            Polynomial::new([0, 0, 1, 0, 0, 0, 0]), // x^2
            Polynomial::new([1, 0, 1, 0, 0, 0, 0]), // x^2 + 1
            Polynomial::new([1, 0, 0, 0, 1, 0, 0]), // x^4 + 1
            Polynomial::new([1, 1, 1, 1, 1, 0, 0]), // x^4 + x^3 + x^2 + x + 1
            Polynomial::new([0, 0, 0, 0, 0, 0, 1]), // x^6
            Polynomial::new([1, 0, 0, 0, 0, 0, 1]), // x^6 + 1
            Polynomial::new([1, 0, 0, 1, 0, 0, 1]), // x^6 + x^3 + 1
            Polynomial::new([2, 5, 9, 3, 7, 0, 1]), // x^6 + 7x^4 + 3x^3 + 9x^2 + 5x + 2
        ] {
            let p_inv = p.inv(13).unwrap();
            assert_eq!(p_inv.mul(&p, 13), one.clone());
        }

        for p in [
            Polynomial::new([1, 0, 0, 1, 0, 0, 0]), // x^3 + 1
            Polynomial::new([1, 0, 0, 0, 0, 1, 0]), // x^5 + 1
        ] {
            let p_inv = p.inv(13);
            assert!(p_inv.is_none());
        }
    }
}
