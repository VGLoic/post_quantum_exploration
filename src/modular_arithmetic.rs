/// Computes a mod n, ensuring the result is non-negative
pub fn modulo(a: i128, n: u64) -> u64 {
    let n_as_i128: i128 = n as i128;
    let mut res = a % n_as_i128;
    if res < 0 {
        res += n_as_i128;
    }
    res as u64
}

/// Computes a + b (mod n)
pub fn modulo_add(a: u64, b: u64, n: u64) -> u64 {
    ((a as u128 + b as u128) % (n as u128)) as u64
}

/// Computes -a (mod n)
pub fn modulo_neg(a: u64, n: u64) -> u64 {
    let a_mod = a % n;
    if a_mod != 0 { n - a_mod } else { 0 }
}

/// Computes a * b (mod n)
pub fn modulo_mul(a: u64, b: u64, n: u64) -> u64 {
    let product = (a as u128) * (b as u128);
    (product % (n as u128)) as u64
}

/// Computes a^(-1) (mod n) using the Extended Euclidean Algorithm
/// Returns None if a has no inverse mod n (i.e. if gcd(a, n) != 1)
pub fn modulo_inv(a: u64, n: u64) -> Option<u64> {
    if a == 0 {
        return None;
    }
    if a == 1 {
        return Some(1);
    }

    let (mut new_r, mut r) = ((a as i128), (n as i128));
    let (mut new_t, mut t) = (1_i128, 0_i128);

    while new_r != 0 {
        let q = r / new_r;
        (new_r, r) = (r - q * new_r, new_r);
        (new_t, t) = (t - q * new_t, new_t);
    }

    if r != 1 {
        return None;
    }

    Some(modulo(t, n))
}

pub fn modulo_exp(a: u64, e: u64, n: u64) -> u64 {
    let mut result = 1;
    let mut base = a;
    let mut exponent = e;
    while exponent > 0 {
        if exponent % 2 == 1 {
            result = modulo_mul(result, base, n);
        }
        exponent >>= 1;
        base = modulo_mul(base, base, n);
    }

    result
}

#[allow(dead_code)]
/// Return greatest common divisor of two integers
///
/// The Euclid's algorithm is used:
/// Let us consider a and b, with a > b
/// Let us consider d, a divisor of a and b, then a = d * alpha, b = d * beta => d also divides a - b = d (alpha - beta).
/// Therefore gcd(a, b) = gcd(a - b, b).
/// We can actually keep decreasing b as long as the first input is greater than b, this is the same as taking the rest of a divided by b, so finally gcd(a, b) = gcd(a % b, b)
///
/// The algorithm goes then as:
///     - identify which input is the larger and which is the smaller,
///     - computes r the rest of the division of larger by smaller,
///         - if zero, the smaller is the gcd,
///         - else, return gcd(smaller, r)
fn gcd(a: u64, b: u64) -> u64 {
    let (mut larger, mut smaller) = if a > b { (a, b) } else { (b, a) };
    if smaller == 0 {
        return larger;
    }
    let mut r = larger % smaller;
    while r != 0 {
        larger = smaller;
        smaller = r;
        r = larger % smaller;
    }
    smaller
}

#[cfg(test)]
mod modular_arithmetic_tests {
    use super::*;

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(0, 1555), 1555);
        assert_eq!(gcd(1, 1555), 1);
        assert_eq!(gcd(132, 0), 132);
        assert_eq!(gcd(132, 1), 1);
        assert_eq!(gcd(27, 33), 3);
        assert_eq!(gcd(1025, 2050), 1025);
        assert_eq!(gcd(2050, 2055), 5);
        assert_eq!(gcd(37, 2), 1);
    }

    #[test]
    fn test_modulo() {
        assert_eq!(modulo(5, 3), 2);
        assert_eq!(modulo(-1, 3), 2);
        assert_eq!(modulo(-4, 3), 2);
        assert_eq!(modulo(10, 7), 3);
        assert_eq!(modulo(-10, 7), 4);
    }

    #[test]
    fn test_modulo_add() {
        assert_eq!(modulo_add(5, 7, 10), 2);
        assert_eq!(modulo_add(0, 0, 13), 0);
        assert_eq!(modulo_add(12, 1, 13), 0);
    }

    #[test]
    fn test_modulo_neg() {
        assert_eq!(modulo_neg(3, 10), 7);
        assert_eq!(modulo_neg(0, 10), 0);
        assert_eq!(modulo_neg(10, 10), 0);
        assert_eq!(modulo_neg(1, 13), 12);
    }

    #[test]
    fn test_modulo_mul() {
        assert_eq!(modulo_mul(3, 4, 5), 2);
        assert_eq!(modulo_mul(0, 10, 7), 0);
        assert_eq!(modulo_mul(6, 6, 7), 1);
        assert_eq!(modulo_mul(5, 5, 5), 0);
        assert_eq!(modulo_mul(1_000_000, 1_000_000, 1_000_000_007), 999993007);
    }

    #[test]
    fn test_modulo_inv_basic() {
        // 3 * 7 = 21 ≡ 1 mod 10, so inv(3, 10) = 7
        assert_eq!(modulo_inv(3, 10), Some(7));
        // 2 has no inverse mod 4
        assert_eq!(modulo_inv(2, 4), None);
        // 1 is always its own inverse
        assert_eq!(modulo_inv(1, 13), Some(1));
        // 0 has no inverse
        assert_eq!(modulo_inv(0, 13), None);
    }

    #[test]
    fn test_modulo_inv_large_prime() {
        let n = 1_000_000_007;
        let a = rand::random();
        let inv = modulo_inv(a, n).unwrap();
        assert_eq!(modulo_mul(a, inv, n), 1);
    }

    #[test]
    fn test_modulo_exp_basic() {
        // 2^3 mod 5 = 8 mod 5 = 3
        assert_eq!(modulo_exp(2, 3, 5), 3);
        // 5^0 mod 7 = 1
        assert_eq!(modulo_exp(5, 0, 7), 1);
        // 0^5 mod 7 = 0
        assert_eq!(modulo_exp(0, 5, 7), 0);
        // 7^1 mod 13 = 7
        assert_eq!(modulo_exp(7, 1, 13), 7);
        // 3^4 mod 5 = 81 mod 5 = 1
        assert_eq!(modulo_exp(3, 4, 5), 1);
    }

    #[test]
    fn test_modulo_exp_large_exponent() {
        // 2^30 mod 1009 = 1073741824 = 1064164 * 1009 + 348
        assert_eq!(modulo_exp(2, 30, 1009), 348);
    }
}
