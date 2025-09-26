#[derive(Debug, PartialEq, Clone, Copy)]
pub struct PrimeFieldElement<const MODULUS: u64>(u64);

impl<const MODULUS: u64> From<u64> for PrimeFieldElement<MODULUS> {
    fn from(value: u64) -> Self {
        PrimeFieldElement::new(value)
    }
}

impl<const MODULUS: u64> PrimeFieldElement<MODULUS> {
    pub fn new(value: u64) -> Self {
        Self(value % MODULUS)
    }

    pub fn inner(&self) -> u64 {
        self.0
    }

    pub fn add(&self, other: &Self) -> Self {
        Self((self.0 + other.0) % MODULUS)
    }

    pub fn mul(&self, other: &Self) -> Self {
        let product = (self.0 as u128) * (other.0 as u128);
        Self((product % MODULUS as u128) as u64)
    }

    pub fn neg(&self) -> Self {
        if self.0 == 0 {
            Self(0)
        } else {
            Self(MODULUS - self.0)
        }
    }

    fn modulo(a: i128) -> Self {
        let mut res = a;
        while res < 0 {
            res += Into::<i128>::into(MODULUS);
        }
        Self(res as u64)
    }

    pub fn inv(&self) -> Option<Self> {
        if self.0 == 0 {
            return None;
        }
        if self.0 == 1 {
            return Some(Self::new(1));
        }

        // At each iteration, we compute the quotient and rest of larger by smaller, i.e. larger = quotient * smaller + rest
        // if rest is not zero, we record the quotient and we repeat with `(smaller, larger % smaller)`
        let mut quotients = vec![];

        let mut larger = MODULUS;
        let mut smaller = self.0;
        let mut r = larger % smaller;

        while r != 0 {
            // The current quotient is registered
            quotients.push(larger / smaller);
            // We prepare the next loop iteration with (smaller, larger % smaller)
            larger = smaller;
            smaller = r;
            // We prepare the next loop iteration by computing the loop condition
            r = larger % smaller;
        }

        // If `smaller != 1`, `gcd(a, N) != 1` and `a` does not have a unique inverse
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

        Some(Self::modulo(beta))
    }
}

impl<const MODULUS: u64> std::fmt::Display for PrimeFieldElement<MODULUS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod prime_field_test {
    use super::*;
    #[test]
    fn test_inverse() {
        type F1B7 = PrimeFieldElement<1_000_000_007>;
        assert!(F1B7::new(0).inv().is_none());
        assert_eq!(F1B7::new(1).inv().unwrap(), F1B7::new(1));
        let p = F1B7::new(rand::random());
        assert_eq!(p.inv().unwrap().mul(&p), F1B7::new(1));
    }
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
mod modular_arithmetic_test {
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

    type FP4 = PrimeFieldElement<4>;
    type FP5 = PrimeFieldElement<5>;
    type FP7 = PrimeFieldElement<7>;
    type FP10 = PrimeFieldElement<10>;
    type FP13 = PrimeFieldElement<13>;

    #[test]
    fn test_modulo_add() {
        assert_eq!(FP10::new(5).add(&7.into()), 2.into());
        assert_eq!(FP13::new(0).add(&0.into()), 0.into());
        assert_eq!(FP13::new(12).add(&1.into()), 0.into());
    }

    #[test]
    fn test_modulo_mul() {
        assert_eq!(FP5::new(3).mul(&4.into()), 2.into());
        assert_eq!(FP7::new(0).mul(&10.into()), 0.into());
        assert_eq!(FP7::new(6).mul(&6.into()), 1.into());
    }

    #[test]
    fn test_modulo_inv_basic() {
        // 3 * 7 = 21 ≡ 1 mod 10, so inv(3, 10) = 7
        assert_eq!(FP10::new(3).inv(), Some(7.into()));
        // 2 has no inverse mod 4
        assert_eq!(FP4::new(2).inv(), None);
        // 1 is always its own inverse
        assert_eq!(FP13::new(1).inv(), Some(1.into()));
        // 0 has no inverse
        assert_eq!(FP13::new(0).inv(), None);
    }

    #[test]
    fn test_modulo_inv_large_prime() {
        type FP1B7 = PrimeFieldElement<1_000_000_007>;
        let a = rand::random();
        let inv = FP1B7::new(a).inv().unwrap();
        assert_eq!(inv.mul(&a.into()), 1.into());
    }
}
