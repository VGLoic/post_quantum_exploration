pub fn modulo(a: i128, n: u64) -> u64 {
    let mut res = a;
    while res < 0 {
        res += Into::<i128>::into(n);
    }
    (res as u64) % n
}

pub fn modulo_add(a: u64, b: u64, n: u64) -> u64 {
    (a + b) % n
}

pub fn modulo_neg(a: u64, n: u64) -> u64 {
    let a = a % n;
    if a == 0 { 0 } else { n - a }
}

pub fn modulo_mul(a: u64, b: u64, n: u64) -> u64 {
    let product = (a as u128) * (b as u128);
    (product % n as u128) as u64
}

pub fn modulo_inv(a: u64, n: u64) -> Option<u64> {
    if a == 0 {
        return None;
    }
    if a == 1 {
        return Some(1);
    }

    // At each iteration, we compute the quotient and rest of larger by smaller, i.e. larger = quotient * smaller + rest
    // if rest is not zero, we record the quotient and we repeat with `(smaller, larger % smaller)`
    let mut quotients = vec![];

    let mut larger = n;
    let mut smaller = modulo(a.into(), n);
    let mut r = larger % smaller;

    while r != 0 {
        // The current quotient is registered
        quotients.push(larger / smaller);

        // We repeat with (smaller, larger % smaller)
        larger = smaller;
        smaller = r;
        r = larger % smaller;
    }

    // If the gcd is not 1, then a has no inverse mod n
    if smaller != 1 {
        return None;
    }

    // At this point we have `smaller` to 1 without rest,
    // going to the last record we then have something of the form larger_k = q_k * smaller_k + 1,
    // we isolate 1 as 1 = larger_k - q_k * smaller_k.
    // Let us define alpha_k and beta_k with the rewrite of the equation: `1 = alpha_k * larger_k + beta_k * smaller_k`,
    // Using the relations: `larger_k = smaller_(k-1)` and `smaller_k = larger_(k-1) - q_(k-1) * smaller_(k-1)`, we can define the suite for `alpha_k` and `beta_k` as:
    // ```
    // alpha_(k-1) = beta_k, alpha_k = 1,
    // beta_(k-1) = alpha_k - beta_k * q_(k-1), beta_k = -q_k
    // ```
    // Now that we have our definitions, we can use the registered quotients in order to recompute the `beta_0`

    // Note that quotients list is necessarily non empty, otherwise `smaller` would be a divider of the prime `N`, therefore `smaller = 1` but this case is treated at the start
    let last_quotient = quotients
        .pop()
        .unwrap_or_else(|| unreachable!("quotient is necessarily non empty"));

    let mut alpha: i128 = 1;
    let mut beta: i128 = -Into::<i128>::into(last_quotient);

    for quotient in quotients.iter().rev() {
        (alpha, beta) = (beta, alpha - beta * Into::<i128>::into(*quotient));
    }

    Some(modulo(beta, n))
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
}
