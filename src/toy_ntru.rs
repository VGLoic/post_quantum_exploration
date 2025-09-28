#[cfg(test)]
mod ntru_tests {
    use crate::modular_arithmetic::*;

    const RING_DEGREE: usize = 7;

    #[derive(Debug, Clone, PartialEq)]
    struct Polynomial {
        // coefficients are in ascending degrees, i.e. a0 + a_1 * x + a_2 * x^2 + ... + a_6 * x^6
        coefficients: [u64; RING_DEGREE],
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

            // let resulting_polynomial_degree = self.degree() - other.degree();

            // x^4 + 2x^2 + 8x + 1 / x^2 = x^2 * (x^2 + 2) + 8x + 1
            // x^4 + 2x^2 + 8x + 1 / (x^2 + 3) = x^2 * (x^2 + 3) - x^2 + 8x + 1 = x^2 * (x^2 + 3) - 1 * (x^2 + 3) + 8x + 4 = (x^2 - 1) * (x^2 + 3) + 8x + 4

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
    }

    // fn gcd(a: &Polynomial, b: &Polynomial) -> Polynomial {
    // }

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
        // x^6 + 2x + 1 = x^3 * (x^3 + 1) -x^3 + 2x + 1 = (x^3 - 1) * (x^3 + 1) + 2x +2
        let (quotient, remainder) = Polynomial::new(p1).div(&Polynomial::new(p2), 13).unwrap();
        assert_eq!(quotient, Polynomial::new([12, 0, 0, 1, 0, 0, 0])); // x^3 - 1
        assert_eq!(remainder, Polynomial::new([2, 2, 0, 0, 0, 0, 0])); // 2x + 2
    }
}
