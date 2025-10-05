use crate::{modular_arithmetic::*, primefield::PrimeFieldElement};

#[derive(Debug, PartialEq, Clone)]
pub struct Polynomial<const N: u32> {
    // Coefficients in ascending order
    coefficients: Vec<PrimeFieldElement<N>>,
}

impl<const N: u32> Default for Polynomial<N> {
    fn default() -> Self {
        Self::new(vec![])
    }
}

impl<const N: u32> Polynomial<N> {
    pub fn new(coefficients: Vec<PrimeFieldElement<N>>) -> Self {
        let mut filtered_coefficients = coefficients;
        let zero = PrimeFieldElement::from(0);
        while filtered_coefficients.last().is_some_and(|&c| c == zero) {
            filtered_coefficients.pop();
        }
        Self {
            coefficients: filtered_coefficients,
        }
    }

    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            return 0;
        }
        self.coefficients.len() - 1
    }

    pub fn add(&self, other: &Self) -> Self {
        let res_degree = self.degree().max(other.degree());
        let mut coefficients = Vec::with_capacity(1 + res_degree);
        let zero = PrimeFieldElement::from(0);
        for i in 0..=res_degree {
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = other.coefficients.get(i).unwrap_or(&zero);
            coefficients.push(a.add(b));
        }
        Self::new(coefficients)
    }

    pub fn neg(&self) -> Self {
        let mut coefficients = Vec::with_capacity(1 + self.degree());
        for c in self.coefficients.iter() {
            coefficients.push(c.neg());
        }
        Self::new(coefficients)
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    pub fn mul(&self, other: &Self) -> Self {
        let res_degree = self.degree() + other.degree();
        let mut coefficients: Vec<PrimeFieldElement<N>> = Vec::with_capacity(1 + res_degree);

        for (i, a) in self.coefficients.iter().enumerate() {
            for (j, b) in other.coefficients.iter().enumerate() {
                let contribution = a.mul(b);
                if let Some(c) = coefficients.get_mut(i + j) {
                    *c = c.add(&contribution)
                } else {
                    coefficients.push(contribution);
                }
            }
        }

        Self::new(coefficients)
    }

    pub fn scalar_mul(&self, a: PrimeFieldElement<N>) -> Self {
        let self_degree = self.degree();
        let mut scaled_coefficients = Vec::with_capacity(self_degree + 1);
        for c in self.coefficients.iter() {
            scaled_coefficients.push(c.mul(&a));
        }
        Self::new(scaled_coefficients)
    }

    pub fn scale_degree(&self, a: usize) -> Self {
        if self.is_zero() {
            return Self::default();
        }
        let res_degree = self.degree() + a;
        let mut res_coefficients = Vec::with_capacity(1 + res_degree);
        for _ in 0..a {
            res_coefficients.push(0.into());
        }
        for i in 0..=self.degree() {
            res_coefficients.push(self.coefficients[i]);
        }
        Self::new(res_coefficients)
    }

    pub fn div(&self, other: &Self) -> Option<(Self, Self)> {
        if other.is_zero() {
            return None;
        }

        let num_degree = self.degree();
        let divider_degree = other.degree();
        if divider_degree > num_degree {
            return Some((Polynomial::default(), self.clone()));
        }

        let quotient_degree = num_degree - divider_degree;
        let mut quotient_coefficients = Vec::with_capacity(1 + quotient_degree);
        let mut remainder = self.clone();

        let divider_leading_coefficient_inv = other.coefficients[divider_degree].inv()?;
        let zero = 0.into();
        for i in (divider_degree..=num_degree).rev() {
            if let Some(remainder_coeff) = remainder.coefficients.get(i) {
                let q = remainder_coeff.mul(&divider_leading_coefficient_inv);
                let to_be_sub = other.scalar_mul(q).scale_degree(i - divider_degree);

                quotient_coefficients.push(q);
                remainder = remainder.add(&to_be_sub.neg());
            } else {
                quotient_coefficients.push(zero);
            }
        }

        remainder = Self::new(remainder.coefficients);
        quotient_coefficients.reverse();

        Some((Self::new(quotient_coefficients), remainder))
    }

    pub fn evaluate(&self, x: PrimeFieldElement<N>) -> PrimeFieldElement<N> {
        let mut evaluation = self
            .coefficients
            .get(0)
            .map(|v| v.clone())
            .unwrap_or(PrimeFieldElement::from(0));
        let mut power_of_x_evaluated = PrimeFieldElement::from(1);
        for c in self.coefficients.iter().skip(1) {
            power_of_x_evaluated = power_of_x_evaluated.mul(&x);
            evaluation = evaluation.add(&c.mul(&power_of_x_evaluated));
        }
        evaluation
    }
}

impl<const N: u32> std::fmt::Display for Polynomial<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let displayed =
            self.coefficients
                .iter()
                .enumerate()
                .fold("".to_string(), |mut state, (i, &v)| {
                    if v != 0.into() {
                        let addendum = match i {
                            0 => v.to_string(),
                            1 => format!("{v}x"),
                            other => format!("{v}x^{other}"),
                        };
                        state += addendum.as_str();
                    }
                    state
                });
        write!(f, "{displayed}")
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    type Polynomial1B7 = Polynomial<1_000_000_007>;

    #[test]
    fn test_creation() {
        let created = Polynomial1B7::new(vec![1.into(), 2.into(), 0.into()]);
        assert_eq!(created, Polynomial1B7::new(vec![1.into(), 2.into()]));
        assert_eq!(created.degree(), 1);

        let created = Polynomial1B7::new(vec![0.into(), 0.into(), 0.into()]);
        assert_eq!(created, Polynomial1B7::new(vec![]));
        assert_eq!(created.degree(), 0);
    }

    #[test]
    fn test_addition() {
        let a = Polynomial1B7::new(vec![1.into(), 3.into(), (1_000_000_007 - 1).into()]);
        let b = Polynomial1B7::new(vec![
            1_000_000_007.into(),
            2.into(),
            2.into(),
            1_000_000_000.into(),
        ]);

        assert_eq!(
            a.add(&b),
            Polynomial1B7::new(vec![1.into(), 5.into(), 1.into(), 1_000_000_000.into()])
        );
    }

    #[test]
    fn test_neg() {
        let mut a: Vec<PrimeFieldElement<1_000_000_007>> = Vec::with_capacity(100);
        for _ in 0..100 {
            let c: u32 = rand::random();
            a.push(c.into())
        }

        let p = Polynomial::new(a);
        assert_eq!(p.add(&p.neg()), Polynomial::new(vec![]));
    }

    #[test]
    fn test_multiplication() {
        let p1 = Polynomial1B7::new([1, 1, 1, 1, 1, 1, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        let p2 = Polynomial1B7::new(
            [2, 4, 8, 16, 32, 64, 128]
                .map(PrimeFieldElement::from)
                .to_vec(),
        ); // 128x^6 + 64x^5 + 32x^4 + 16x^3 + 8x^2 + 4x + 2
        let expected = Polynomial1B7::new(
            [2, 6, 14, 30, 62, 126, 254, 252, 248, 240, 224, 192, 128]
                .map(PrimeFieldElement::from)
                .to_vec(),
        );
        assert_eq!(p1.mul(&p2), expected);
    }

    #[test]
    fn test_evaluation() {
        let p = Polynomial1B7::new([1, 1, 1, 1, 1, 1, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        assert_eq!(p.evaluate(1.into()), 7.into());
        assert_eq!(p.evaluate(0.into()), 1.into());
        let p = Polynomial1B7::new([2, 1_000_000_005, 4].map(PrimeFieldElement::from).to_vec());
        assert_eq!(p.evaluate(0.into()), 2.into());
        assert_eq!(p.evaluate(1.into()), 4.into());
        assert_eq!(p.evaluate(4.into()), 58.into()); 
    }

    #[test]
    fn test_division() {
        let p1 = Polynomial1B7::new([0, 0, 0, 1, 0, 0, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + x^3
        let p2 = Polynomial1B7::new([1, 0, 0, 1].map(PrimeFieldElement::from).to_vec()); // x^3 + 1
        let (quotient, remainder) = p1.div(&p2).unwrap();
        assert_eq!(
            quotient,
            Polynomial1B7::new([0, 0, 0, 1].map(PrimeFieldElement::from).to_vec())
        ); // x^3
        assert_eq!(remainder, Polynomial1B7::new(vec![])); // 0

        let p1 = Polynomial1B7::new([1, 2, 0, 0, 0, 0, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + 2x + 1
        let p2 = Polynomial1B7::new([1, 0, 0, 1].map(PrimeFieldElement::from).to_vec()); // x^3 + 1
        // x^6 + 2x + 1 = x^3 * (x^3 + 1) -x^3 + 2x + 1 = (x^3 - 1) * (x^3 + 1) + 2x + 2
        let (quotient, remainder) = p1.div(&p2).unwrap();
        assert_eq!(
            quotient,
            Polynomial1B7::new(
                [-1, 0, 0, 1]
                    .map(|v| PrimeFieldElement::from(modulo(v, 1_000_000_007)))
                    .to_vec()
            )
        ); // x^3 - 1
        assert_eq!(remainder, Polynomial1B7::new(vec![2.into(), 2.into()])); // 2x + 2
    }
}
