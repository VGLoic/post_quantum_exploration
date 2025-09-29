#![allow(dead_code)]

use super::modular_arithmetic;

/// Computes a + b (mod n) in centered coordinates, i.e. between -n/2 and n/2
pub fn modulo_add(a: i64, b: i64, n: u32) -> i64 {
    to_centered_coordinates(Into::<i128>::into(a) + Into::<i128>::into(b), n)
}

/// Computes -a (mod n) in centered coordinates, i.e. between -n/2 and n/2
pub fn modulo_neg(a: i64, n: u32) -> i64 {
    -to_centered_coordinates(a.into(), n)
}

/// Computes a * b (mod n) in centered coordinates, i.e. between -n/2 and n/2
pub fn modulo_mul(a: i64, b: i64, n: u32) -> i64 {
    let product = Into::<i128>::into(a) * Into::<i128>::into(b);
    to_centered_coordinates(product, n)
}

/// Computes a^(-1) (mod n) in centered coordinates, i.e. between -n/2 and n/2
pub fn modulo_inv(a: i64, n: u32) -> Option<i64> {
    modular_arithmetic::modulo_inv(modular_arithmetic::modulo(a.into(), n), n)
        .map(|inv| to_centered_coordinates(inv.into(), n))
}

pub fn to_centered_coordinates(a: i128, n: u32) -> i64 {
    let res = (a % (Into::<i128>::into(n))) as i64;
    let half_n = Into::<i64>::into(n / 2);
    if res < -half_n {
        res + (Into::<i64>::into(n))
    } else if res > half_n {
        res - (Into::<i64>::into(n))
    } else {
        res
    }
}

#[cfg(test)]
mod modular_arithmetic_tests {
    use super::*;

    #[test]
    fn test_modulo() {
        assert_eq!(to_centered_coordinates(5, 3), -1);
        assert_eq!(to_centered_coordinates(4, 3), 1);
        assert_eq!(to_centered_coordinates(3, 3), 0);
        assert_eq!(to_centered_coordinates(-1, 3), -1);
        assert_eq!(to_centered_coordinates(-4, 3), -1);
        assert_eq!(to_centered_coordinates(10, 7), 3);
        assert_eq!(to_centered_coordinates(-10, 7), -3);
    }

    #[test]
    fn test_modulo_add() {
        assert_eq!(modulo_add(5, 7, 10), 2);
        assert_eq!(modulo_add(6, 9, 10), 5);
        assert_eq!(modulo_add(8, 9, 10), -3);
        assert_eq!(modulo_add(0, 0, 13), 0);
        assert_eq!(modulo_add(12, 1, 13), 0);
    }

    #[test]
    fn test_modulo_neg() {
        assert_eq!(modulo_neg(3, 10), -3);
        assert_eq!(modulo_neg(0, 10), 0);
        assert_eq!(modulo_neg(10, 10), 0);
        assert_eq!(modulo_neg(1, 13), -1);
        assert_eq!(modulo_neg(9, 13), 4);
    }

    #[test]
    fn test_modulo_mul() {
        assert_eq!(modulo_mul(3, 4, 5), 2);
        assert_eq!(modulo_mul(0, 10, 7), 0);
        assert_eq!(modulo_mul(6, 6, 7), 1);
        assert_eq!(modulo_mul(5, 5, 5), 0);
        assert_eq!(modulo_mul(4, 8, 3), -1);
    }

    #[test]
    fn test_modulo_inv_basic() {
        // 3 * 7 = 21 ≡ 1 mod 10, so inv(3, 10) = 7 => -3 in centered
        assert_eq!(modulo_inv(3, 10), Some(-3));
        // 2 has no inverse mod 4
        assert_eq!(modulo_inv(2, 4), None);
        // 1 is always its own inverse
        assert_eq!(modulo_inv(1, 13), Some(1));
        // 0 has no inverse
        assert_eq!(modulo_inv(0, 13), None);
    }

    #[test]
    fn test_modulo_inv_prime() {
        let n = 19;
        let mut a = rand::random();
        if a % (n as i64) == 0 {
            a -= 1;
        }
        let inv = modulo_inv(a, n).unwrap();
        assert_eq!(modulo_mul(a, inv, n), 1);
    }
}
