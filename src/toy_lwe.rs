#[cfg(test)]
mod toy_lwe_tests {
    use crate::primefield::*;

    /*
     * This is a toy example for Learning with Errors (LWE), not optimized for performance or size.
     * It illustrates the basic principles of LWE-based key exchange using simple vector and matrix operations.
     *
     * The code defines a small 4x4 matrix and 4-dimensional vectors, simulating the LWE key exchange process.
     * An assumption is made that the matrix is symmetric to simplify the example. In practice, LWE implementations use much larger dimensions and more complex error distributions.
     *
     * The key exchange process involves two parties, Alice and Bob, who each generate a secret vector and an error vector with low amplitude.
     * They compute their public keys by multiplying the shared matrix with their secret vectors and adding their respective error vectors.
     * Both parties then derive a shared secret by computing the scalar product of their secret vector with the other party's public key, followed by truncation to reduce noise.
     * The truncation step is crucial in LWE to ensure that both parties arrive at the same shared secret despite the presence of noise.
     */

    type FieldElement = PrimeFieldElement<1_000_000_007>;

    #[derive(Debug, PartialEq)]
    struct Vector {
        v: [FieldElement; 4],
    }

    impl Vector {
        fn new(v: [FieldElement; 4]) -> Self {
            Self { v }
        }

        fn scalar_product(&self, b: &Vector) -> FieldElement {
            self.v
                .iter()
                .zip(&b.v)
                .fold(FieldElement::new(0), |acc, pair| {
                    acc.add(&pair.0.mul(pair.1))
                })
        }

        fn add(&self, b: &Vector) -> Vector {
            let mut result = [FieldElement::new(0); 4];
            for (i, pair) in self.v.iter().zip(&b.v).enumerate() {
                result[i] = pair.0.add(pair.1);
            }
            Vector::new(result)
        }
    }

    impl std::fmt::Display for Vector {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "[{}]", &self.v.map(|v| v.to_string()).join(", "))
        }
    }

    struct Matrix {
        rows: [[FieldElement; 4]; 4],
    }

    impl Matrix {
        fn new(rows: [[FieldElement; 4]; 4]) -> Self {
            Self { rows }
        }

        fn mult(&self, v: &Vector) -> Vector {
            let mut result = [FieldElement::new(0); 4];
            for (i, row) in self.rows.iter().enumerate() {
                result[i] = row
                    .iter()
                    .zip(&v.v)
                    .fold(FieldElement::new(0), |acc, pair| {
                        acc.add(&pair.0.mul(pair.1))
                    });
            }
            Vector::new(result)
        }
    }

    impl std::fmt::Display for Matrix {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let mut displayed = "".to_string();
            for row in &self.rows {
                displayed += format!("[{}]\n", row.map(|v| v.to_string()).join(", ")).as_str();
            }
            write!(f, "\n{displayed}\n")
        }
    }

    fn truncate(value: FieldElement, order: u8) -> FieldElement {
        let truncate_factor = 10u32.pow(order.into());
        let truncated_value = value.inner() / truncate_factor * truncate_factor;
        FieldElement::new(truncated_value)
    }

    #[test]
    fn test_vector_ops() {
        let one = FieldElement::new(1);
        let zero = FieldElement::new(0);
        let a = Vector::new([one, zero, zero, zero]);
        let b = Vector::new([zero, one, zero, zero]);
        assert_eq!(a.add(&b), Vector::new([one, one, zero, zero]));
        assert_eq!(a.scalar_product(&b), zero);
        let m = Matrix::new([
            [one, zero, zero, zero],
            [zero, one, zero, zero],
            [zero, zero, one, zero],
            [zero, zero, zero, one],
        ]);
        assert_eq!(m.mult(&a), a);
        assert_eq!(m.mult(&b), b);
    }

    fn low_amplitude_random_vector(order: u8) -> Vector {
        let mut v = [FieldElement::new(0); 4];
        let higher_limit = 10u32.pow(order.into());
        for value in v.iter_mut() {
            *value = FieldElement::new(rand::random_range(0..higher_limit));
        }
        Vector { v }
    }

    #[test]
    fn test_lwe_dh() {
        let mut a = Matrix::new(
            (0..4)
                .map(|_| {
                    let mut row = [FieldElement::new(0); 4];
                    for item in row.iter_mut() {
                        *item = FieldElement::new(rand::random_range(1_000_000_000_u32..u32::MAX));
                    }
                    row
                })
                .collect::<Vec<[FieldElement; 4]>>()
                .try_into()
                .unwrap(),
        );
        // A is taken as symmetric
        for row_index in 1..4 {
            for column_index in 0..row_index {
                a.rows[row_index][column_index] = a.rows[column_index][row_index];
            }
        }

        let vector_amplitude = 2;

        let s_alice = low_amplitude_random_vector(vector_amplitude);
        let e_alice = low_amplitude_random_vector(vector_amplitude);
        // p_a = A * s_a + e_a
        let public_key_alice = a.mult(&s_alice).add(&e_alice);

        let s_bob = low_amplitude_random_vector(vector_amplitude);
        let e_bob = low_amplitude_random_vector(vector_amplitude);
        // p_b = A * s_b + e_b
        let public_key_bob = a.mult(&s_bob).add(&e_bob);

        // This should be just enough to cover the scalar product of two low amplitude vectors
        let truncation_order = 1 + 1 + vector_amplitude * vector_amplitude;

        // s_b . p_b = s_b . (A * s_a + e_a) = s_b . (A * s_a) + s_b . e_a = s_b^transposed * A * s_a + s_b . e_a
        let k_bob = truncate(public_key_alice.scalar_product(&s_bob), truncation_order);
        // s_a . p_b = s_a . (A * s_b + e_b) = s_a . (A * s_b) + s_a . e_b = s_a^transposed * A * s_b + s_a . e_b = s_b^transposed * A * s_a + s_a * e_b
        let k_alice = truncate(public_key_bob.scalar_product(&s_alice), truncation_order);

        assert_eq!(k_bob, k_alice);
    }
}
