use crate::modular_arithmetic::*;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct PrimeFieldElement<const MODULUS: u32>(u32);

impl<const MODULUS: u32> From<u32> for PrimeFieldElement<MODULUS> {
    fn from(value: u32) -> Self {
        PrimeFieldElement::new(value)
    }
}

impl<const MODULUS: u32> PrimeFieldElement<MODULUS> {
    pub fn new(value: u32) -> Self {
        Self(modulo(value.into(), MODULUS))
    }

    pub fn inner(&self) -> u32 {
        self.0
    }

    pub fn add(&self, other: &Self) -> Self {
        Self(modulo_add(self.0, other.0, MODULUS))
    }

    pub fn mul(&self, other: &Self) -> Self {
        Self(modulo_mul(self.0, other.0, MODULUS))
    }

    pub fn neg(&self) -> Self {
        Self(modulo_neg(self.0, MODULUS))
    }

    pub fn inv(&self) -> Option<Self> {
        modulo_inv(self.0, MODULUS).map(Self)
    }

    pub fn exp(&self, e: u32) -> Self {
        Self(modulo_exp(self.0, e, MODULUS))
    }
}

impl<const MODULUS: u32> std::fmt::Display for PrimeFieldElement<MODULUS> {
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

#[cfg(test)]
mod modular_arithmetic_test {
    use super::*;

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
