use crate::primefield::PrimeFieldElement;

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

    pub fn interpolate_from_roots(roots: Vec<PrimeFieldElement<N>>) -> Self {
        Self::interpolate_from_roots_slice(&roots)
    }
    pub fn interpolate_from_roots_slice(roots: &[PrimeFieldElement<N>]) -> Self {
        if roots.is_empty() {
            return Self::default();
        }
        let mut coefficients = Vec::with_capacity(roots.len());
        coefficients.push(PrimeFieldElement::from(1));
        for (i, x) in roots.iter().enumerate() {
            if i % 10_000 == 0 {
                println!("[interpolation from roots] reached a 10_000");
            }
            let x_neg = x.neg();
            // Handle leading coefficient
            coefficients.push(1.into());
            // For each coefficient from 1 to i (incuded) we have coeff[j] = coeff[j - 1] - x * coeff[j]
            for j in (1..=i).rev() {
                coefficients[j] = coefficients[j - 1].add(&coefficients[j].mul(&x_neg));
            }
            // For 0, we have coeff[0] = -x * coeff[0]
            coefficients[0] = coefficients[0].mul(&x_neg);
        }

        Self::new(coefficients)
    }

    pub fn interpolate_from_coordinates(
        points: Vec<PrimeFieldElement<N>>,
        values: Vec<PrimeFieldElement<N>>,
    ) -> Option<Self> {
        if points.is_empty() || points.len() != values.len() {
            return Some(Self::default());
        }

        let master_numerator = Self::interpolate_from_roots_slice(&points);

        let mut p_coefficients = vec![PrimeFieldElement::<N>::from(0); points.len()];
        for (x_i, y_i) in points.into_iter().zip(values) {
            let (numerator, _) = master_numerator.div(&Self::new(vec![x_i.neg(), 1.into()]))?;
            let factor = y_i.mul(&numerator.evaluate(x_i).inv()?);
            for (j, c_j) in numerator.coefficients.into_iter().enumerate() {
                p_coefficients[j] = p_coefficients[j].add(&c_j.mul(&factor));
            }
        }
        Some(Self::new(p_coefficients))
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

    pub fn div(&self, other: &Self) -> Option<(Self, Self)> {
        if other.is_zero() {
            return None;
        }

        let num_degree = self.degree();
        let denominator_degree = other.degree();
        if denominator_degree > num_degree {
            return Some((Polynomial::default(), self.clone()));
        }

        let quotient_degree = num_degree - denominator_degree;
        let mut quotient_coefficients = vec![PrimeFieldElement::from(0); 1 + quotient_degree];
        let mut remainder_coefficients = self.coefficients.clone();

        let divider_leading_coefficient_inv = other.coefficients[denominator_degree].inv()?;
        for i in (denominator_degree..=num_degree).rev() {
            if remainder_coefficients[i].inner() != 0 {
                let scaling_degree = i - denominator_degree;
                let q = remainder_coefficients[i].mul(&divider_leading_coefficient_inv);
                quotient_coefficients[scaling_degree] = q;
                if q.inner() != 0 {
                    let q_neg = q.neg();
                    for (j, denominator_coeff) in other.coefficients.iter().enumerate() {
                        remainder_coefficients[j + scaling_degree] = remainder_coefficients
                            [j + scaling_degree]
                            .add(&denominator_coeff.mul(&q_neg));
                    }
                }
            }
            remainder_coefficients.pop();
        }

        Some((
            Self::new(quotient_coefficients),
            Self::new(remainder_coefficients),
        ))
    }

    pub fn evaluate(&self, x: PrimeFieldElement<N>) -> PrimeFieldElement<N> {
        let mut evaluation = self
            .coefficients
            .first()
            .copied()
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
    use crate::modular_arithmetic::modulo;

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
    fn test_interpolation_from_roots() {
        let number_of_points = rand::random_range(2..=100);
        let points: Vec<PrimeFieldElement<1_000_000_007>> = (0..number_of_points)
            .map(|_| {
                let v: u32 = rand::random();
                PrimeFieldElement::from(v)
            })
            .collect();
        let p = Polynomial1B7::interpolate_from_roots(points.clone());
        for point in points {
            assert_eq!(p.evaluate(point), 0.into());
        }
    }

    #[test]
    fn test_interpolation_from_coordinates() {
        let number_of_points: u32 = rand::random_range(2..=100);
        let points: Vec<PrimeFieldElement<1_000_000_007>> =
            (0..number_of_points).map(PrimeFieldElement::from).collect();
        let values: Vec<PrimeFieldElement<1_000_000_007>> = (0..number_of_points)
            .map(|_| {
                let y: u32 = rand::random();
                PrimeFieldElement::from(y)
            })
            .collect();
        let p =
            Polynomial1B7::interpolate_from_coordinates(points.clone(), values.clone()).unwrap();
        for (x, y) in points.into_iter().zip(values) {
            assert_eq!(p.evaluate(x), y);
        }
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
