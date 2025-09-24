#[cfg(test)]
mod toy_lwe_tests {

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

    #[derive(Debug, PartialEq)]
    struct Vector {
        v: [u128; 4],
    }

    impl Vector {
        fn new(v: [u128; 4]) -> Self {
            Self { v }
        }

        fn scalar_product(&self, b: &Vector) -> u128 {
            self.v.iter().zip(&b.v).fold(0u128, |acc, pair| {
                let contribution = pair.0.checked_mul(*pair.1).unwrap();
                acc.checked_add(contribution).unwrap()
            })
        }

        fn add(&self, b: &Vector) -> Vector {
            let mut result = [0u128; 4];
            for (i, pair) in self.v.iter().zip(&b.v).enumerate() {
                result[i] = pair.0.checked_add(*pair.1).unwrap();
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
        rows: [[u128; 4]; 4],
    }

    impl Matrix {
        fn new(rows: [[u128; 4]; 4]) -> Self {
            Self { rows }
        }

        fn mult(&self, v: &Vector) -> Vector {
            let mut result = [0u128; 4];
            for (i, row) in self.rows.iter().enumerate() {
                result[i] = row.iter().zip(&v.v).fold(0u128, |acc, pair| {
                    let contribution = pair.0.checked_mul(*pair.1).unwrap();
                    acc.checked_add(contribution).unwrap()
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

    fn truncate(value: u128, order: u8) -> u128 {
        let truncate_factor = 10u128.pow(order.into());
        value / truncate_factor * truncate_factor
    }

    #[test]
    fn test_vector_ops() {
        let a = Vector::new([1, 0, 0, 0]);
        let b = Vector::new([0, 1, 0, 0]);
        assert_eq!(a.add(&b), Vector::new([1, 1, 0, 0]));
        assert_eq!(a.scalar_product(&b), 0);
        let m = Matrix::new([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
        assert_eq!(m.mult(&a), a);
        assert_eq!(m.mult(&b), b);
    }

    fn low_amplitude_random_vector(order: u8) -> Vector {
        let mut v = [0u128; 4];
        let higher_limit = 10u128.pow(order.into());
        for value in v.iter_mut() {
            *value = rand::random_range(0..higher_limit)
        }
        Vector { v }
    }

    #[test]
    fn test_lwe_dh() {
        let mut a = Matrix::new(
            (0..4)
                .map(|_| {
                    let mut row = [0u128; 4];
                    for item in row.iter_mut() {
                        *item = rand::random_range(1_000_000_000_u128..1_000_000_000_000_u128);
                    }
                    row
                })
                .collect::<Vec<[u128; 4]>>()
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
