#[cfg(test)]
mod ntru_tests {
    use crate::modular_arithmetic::*;

    const RING_DEGREE: usize = 7;

    #[derive(Debug, Clone, PartialEq, Default)]
    struct Polynomial {
        // coefficients are in ascending degrees, i.e. a0 + a_1 * x + a_2 * x^2 + ... + a_6 * x^6
        coefficients: [u64; RING_DEGREE],
    }

    impl std::fmt::Display for Polynomial {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let terms: Vec<String> = self
                .coefficients
                .iter()
                .enumerate()
                .filter(|(_, coeff)| **coeff != 0)
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
            if terms.is_empty() {
                return write!(f, "0");
            }
            write!(f, "{}", terms.join(" + "))
        }
    }

    impl Polynomial {
        fn new(coefficients: [u64; RING_DEGREE]) -> Self {
            Self { coefficients }
        }

        fn add(&self, other: &Self, n: u64) -> Self {
            let mut result_coefficients = [0u64; RING_DEGREE];

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

        fn neg(&self, n: u64) -> Self {
            let mut coefficients = [0u64; RING_DEGREE];

            for (i, coeff) in self.coefficients.iter().enumerate() {
                coefficients[i] = modulo_neg(*coeff, n);
            }

            Self::new(coefficients)
        }

        fn is_zero(&self) -> bool {
            !self.coefficients.iter().any(|v| *v != 0)
        }

        fn degree(&self) -> usize {
            for i in (0..RING_DEGREE).rev() {
                if self.coefficients[i] != 0 {
                    return i;
                }
            }
            0
        }

        fn mul(&self, other: &Self, n: u64) -> Self {
            let mut raw_coefficients = [0u64; 2 * RING_DEGREE - 1];
            for (i, a_i) in self.coefficients.iter().enumerate() {
                for (j, b_j) in other.coefficients.iter().enumerate() {
                    raw_coefficients[i + j] =
                        modulo_add(raw_coefficients[i + j], modulo_mul(*a_i, *b_j, n), n);
                }
            }

            let (coefficients, to_be_filtered) = raw_coefficients.split_at_mut(RING_DEGREE);

            for (i, coeff) in to_be_filtered.iter().enumerate() {
                // We subtract `coeff * x^i * (x^RING_DEGREE + 1)` to the full polynomial
                // The `RING_DEGREE + i` coefficient will then vanish
                // We only register the subtraction of the `i`th coefficient by `coeff`
                coefficients[i] = modulo_add(coefficients[i], modulo_neg(*coeff, n), n);
            }

            Self::new(coefficients.try_into().unwrap_or_else(|_| {
                unreachable!("coefficients is necessarily of length RING_DEGREE")
            }))
        }

        fn div(&self, other: &Self, n: u64) -> Option<(Self, Self)> {
            let other_degree = other.degree();
            let self_degree = self.degree();
            if self.is_zero() || other.is_zero() || other_degree > self_degree {
                return Some((Polynomial::new([0u64; RING_DEGREE]), self.clone()));
            }

            // Starting from the highest degree term of self, we will subtract multiples of other until we reach the degree of other
            // e.g. if self has degree 6 and other has degree 3, we will do the iterationss for degrees 6, 5, 4 and 3.

            let mut remainder = self.clone();
            let mut quotient_coefficients = [0u64; RING_DEGREE];

            let other_leading_coefficient = other.coefficients[other_degree];

            for i in (other_degree..=self_degree).rev() {
                let other_leading_coefficient_inv = modulo_inv(other_leading_coefficient, n)?;
                let quotient_coefficient =
                    modulo_mul(remainder.coefficients[i], other_leading_coefficient_inv, n);
                quotient_coefficients[i - other_degree] = quotient_coefficient;
                let mut to_be_subtracted_coeff = [0u64; RING_DEGREE];
                to_be_subtracted_coeff[i - other_degree] = quotient_coefficient;
                let to_be_subtracted = Polynomial::new(to_be_subtracted_coeff).mul(other, n).neg(n);
                remainder = remainder.add(&to_be_subtracted, n);
            }

            Some((Polynomial::new(quotient_coefficients), remainder))
        }

        fn inv(&self, n: u64) -> Option<Polynomial> {
            if self.is_zero() {
                return None;
            }
            let self_degree = self.degree();

            if self_degree == 0 {
                let inv_coeff = modulo_inv(self.coefficients[0], n)?;
                return Some(Polynomial::new([inv_coeff, 0, 0, 0, 0, 0, 0]));
            }

            // We can't handle `R(x) = x^RING_DEGREE + 1` in our struct so we will perform a first division step of `R` by `self` in order to remove the x^RING_DEGREE coefficient
            // We end up with larger = R - correction * self
            let (corrected, correction) = {
                let self_leading_coefficient = self.coefficients[self_degree];
                let quotient_coefficient = modulo_inv(self_leading_coefficient, n)?;

                let mut to_be_subtracted_coeff = [0u64; RING_DEGREE];
                to_be_subtracted_coeff[RING_DEGREE - self_degree] = quotient_coefficient;
                let correction = Polynomial::new(to_be_subtracted_coeff);

                // Every coefficient is multiplied by `quotient_coefficient`, increase its order of `RING_DEGREE - self_degree` and is negated, the highest one is ignored
                let degree_increase = RING_DEGREE - self_degree;
                let mut corrected_coefficients = [0u64; RING_DEGREE];
                for (i, coeff) in self.coefficients.iter().enumerate() {
                    if i + degree_increase < RING_DEGREE {
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

    fn gcd(a: &Polynomial, b: &Polynomial, n: u64) -> Option<Polynomial> {
        if a.is_zero() {
            return Some(b.clone());
        }
        if b.is_zero() {
            return Some(a.clone());
        }

        let (larger, mut smaller) = if a.degree() > b.degree() {
            (a, b.clone())
        } else {
            (b, a.clone())
        };
        let (_, mut r) = larger.div(&smaller, n)?;

        // As long as the remainder is not zero, we use the fact that gcd(p, q) = gcd(q, p % q)
        // At some point, either we hit a non zero degree polynomial, either we arrive at order 0 where the remainder will always be zero
        while !r.is_zero() {
            let larger = smaller;
            smaller = r;
            (_, r) = larger.div(&smaller, n)?;
        }

        Some(smaller)
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new([1, 2, 3, 4, 5, 6, 7]);
        let p2 = Polynomial::new([7, 6, 5, 4, 3, 2, 1]);
        let expected = Polynomial::new([8, 8, 8, 8, 8, 8, 8]);
        assert_eq!(p1.add(&p2, 13), expected);

        let p3 = Polynomial::new([12, 12, 12, 12, 12, 12, 12]);
        let p4 = Polynomial::new([1, 1, 1, 1, 1, 1, 1]);
        let expected2 = Polynomial::new([0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(p3.add(&p4, 13), expected2);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new([1, 1, 1, 1, 1, 1, 1]);
        let p2 = Polynomial::new([2, 4, 8, 16, 32, 64, 128]);
        let expected = Polynomial::new([
            999_999_757,
            999_999_765,
            999_999_781,
            999_999_813,
            999_999_877,
            1_000_000_005,
            254,
        ]);
        assert_eq!(p1.mul(&p2, 1_000_000_007), expected);
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
        assert_eq!(quotient, Polynomial::new([12, 0, 0, 1, 0, 0, 0])); // x^3 - 1
        assert_eq!(remainder, Polynomial::new([2, 2, 0, 0, 0, 0, 0])); // 2x + 2
    }

    #[test]
    fn test_polynomial_gcd() {
        let p1 = Polynomial::new([0, 0, 0, 1, 0, 0, 1]); // x^6 + x^3
        let p2 = Polynomial::new([1, 0, 0, 1, 0, 0, 0]); // x^3 + 1
        let g = gcd(&p1, &p2, 13).unwrap();
        assert_eq!(g, p2); // x^3 + 1

        let p1 = Polynomial::new([1, 2, 0, 0, 0, 0, 1]); // x^6 + 2x + 1
        let p2 = Polynomial::new([1, 0, 0, 1, 0, 0, 0]); // x^3 + 1
        let g = gcd(&p1, &p2, 13).unwrap(); // 2(x + 1)
        assert_eq!(g, Polynomial::new([2, 2, 0, 0, 0, 0, 0]));

        let p1 = Polynomial::new([1, 2, 0, 0, 0, 5, 1]); // x^6 + 5x^5 + 2x + 1
        let p2 = Polynomial::new([1, 1, 0, 0, 0, 0, 0]); // x + 1
        let g = gcd(&p1, &p2, 13).unwrap(); // 8
        assert_eq!(g, Polynomial::new([8, 0, 0, 0, 0, 0, 0]));
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
