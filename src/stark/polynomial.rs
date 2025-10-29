use anyhow::anyhow;

use crate::{primefield::PrimeFieldElement, stark::fft};

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

    pub fn interpolate_from_roots(roots: &[PrimeFieldElement<N>]) -> Self {
        if roots.is_empty() {
            return Self::default();
        }
        let mut coefficients = Vec::with_capacity(roots.len());
        coefficients.push(PrimeFieldElement::from(1));
        for (i, x) in roots.iter().enumerate() {
            let x_neg = x.neg();
            // Handle leading coefficient
            coefficients.push(1.into());
            // For each coefficient from 1 to i (included) we have coeff[j] = coeff[j - 1] - x * coeff[j]
            for j in (1..=i).rev() {
                coefficients[j] = coefficients[j - 1].add(&coefficients[j].mul(&x_neg));
            }
            // For 0, we have coeff[0] = -x * coeff[0]
            coefficients[0] = coefficients[0].mul(&x_neg);
        }

        Self::new(coefficients)
    }

    pub fn interpolate_and_evaluate_zpoly(
        roots: impl Iterator<Item = u32>,
        x: &PrimeFieldElement<N>,
    ) -> PrimeFieldElement<N> {
        let mut is_empty = true;
        let mut result = PrimeFieldElement::<N>::from(1);
        for root in roots {
            is_empty = false;
            let neg_root = PrimeFieldElement::<N>::from(N - root);
            result = result.mul(&x.add(&neg_root));
        }
        if is_empty {
            return 0.into();
        }

        result
    }

    pub fn interpolate_and_evaluate_from_coordinates(
        points: &[PrimeFieldElement<N>],
        values: &[PrimeFieldElement<N>],
        x: &PrimeFieldElement<N>,
    ) -> Option<PrimeFieldElement<N>> {
        if points.is_empty() || points.len() != values.len() {
            return Some(0.into());
        }

        let mut result = PrimeFieldElement::<N>::from(0);
        for (i, (x_i, y_i)) in points.iter().zip(values).enumerate() {
            let mut contribution = *y_i;
            for (j, root_j) in points.iter().enumerate() {
                if i != j {
                    let root_j_neg = PrimeFieldElement::new(N - root_j.inner());
                    let numerator = x.add(&root_j_neg);
                    let denominator = x_i.add(&root_j_neg);
                    contribution = contribution.mul(&numerator.mul(&denominator.inv()?));
                }
            }
            result = result.add(&contribution);
        }

        Some(result)
    }

    pub fn interpolate_from_coordinates(
        points: &[PrimeFieldElement<N>],
        values: &[PrimeFieldElement<N>],
    ) -> Result<Self, anyhow::Error> {
        if points.is_empty() || points.len() != values.len() {
            return Ok(Self::default());
        }

        let master_numerator = Self::interpolate_from_roots(points);

        let mut p_coefficients = vec![PrimeFieldElement::<N>::from(0); points.len()];
        for (x_i, y_i) in points.iter().zip(values) {
            let (numerator, _) = master_numerator
                .div(&Self::new(vec![x_i.neg(), 1.into()]))
                .ok_or(anyhow!("unable to divide master numerator with root {x_i}"))?;
            let numerator_ev_inv = numerator.evaluate(x_i).inv().ok_or(anyhow!("failed to invert numerator evaluation; this indicates duplicate points in input, which is not supported"))?;
            let factor = y_i.mul(&numerator_ev_inv);
            for (j, c_j) in numerator.coefficients.into_iter().enumerate() {
                p_coefficients[j] = p_coefficients[j].add(&c_j.mul(&factor));
            }
        }
        Ok(Self::new(p_coefficients))
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

    pub fn mul_by_scalar(&self, a: &PrimeFieldElement<N>) -> Self {
        let mut coefficients: Vec<PrimeFieldElement<N>> = Vec::with_capacity(1 + self.degree());
        for c in &self.coefficients {
            coefficients.push(c.mul(a));
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

    pub fn evaluate(&self, x: &PrimeFieldElement<N>) -> PrimeFieldElement<N> {
        let mut evaluation = self
            .coefficients
            .first()
            .copied()
            .unwrap_or(PrimeFieldElement::from(0));
        let mut power_of_x_evaluated = PrimeFieldElement::from(1);
        for c in self.coefficients.iter().skip(1) {
            power_of_x_evaluated = power_of_x_evaluated.mul(x);
            evaluation = evaluation.add(&c.mul(&power_of_x_evaluated));
        }
        evaluation
    }

    /// Evaluates the polynomial over a whole set of units generated as successive power of a generator.
    ///
    /// The implementation uses the fast Fourier transform, therefore the units set must have a length as 2^k.
    ///
    /// The polynomial must have a degree equal or lower than the length of the units set.
    ///
    /// # Arguments
    /// * `units` - Set of units, generated as successive power of a generator. The set must have a length of 2^k.
    ///
    /// Returns the evaluation of the polynomial over the whole set.
    pub fn fft_evaluate(&self, units: &[PrimeFieldElement<N>]) -> Vec<PrimeFieldElement<N>> {
        fft::fft(&self.coefficients, units)
    }

    pub fn evaluate_as_binomial(
        &self,
        x: &PrimeFieldElement<N>,
        y: &PrimeFieldElement<N>,
        exponent: u32,
    ) -> PrimeFieldElement<N> {
        let mut evaluation = self
            .coefficients
            .first()
            .copied()
            .unwrap_or(PrimeFieldElement::from(0));
        let mut base = PrimeFieldElement::<N>::from(1);
        let mut power_of_x_evaluated = PrimeFieldElement::<N>::from(1);
        for (i, c) in self.coefficients.iter().skip(1).enumerate() {
            if (i + 1) % (exponent as usize) == 0 {
                base = base.mul(y);
                power_of_x_evaluated = base;
            } else {
                power_of_x_evaluated = power_of_x_evaluated.mul(x);
            }
            evaluation = evaluation.add(&c.mul(&power_of_x_evaluated));
        }
        evaluation
    }

    pub fn partially_evaluate_as_binomial(
        &self,
        x: &PrimeFieldElement<N>,
        exponent: usize,
    ) -> Self {
        let mut coefficients = Vec::with_capacity(self.degree() / exponent);

        let mut x_powered = PrimeFieldElement::<N>::from(1);
        for (i, c) in self.coefficients.iter().enumerate() {
            if i % exponent == 0 {
                coefficients.push(*c);
                x_powered = PrimeFieldElement::<N>::from(1)
            } else {
                let index = i / exponent;
                x_powered = x_powered.mul(x);
                coefficients[index] = coefficients[index].add(&c.mul(&x_powered));
            }
        }

        Self::new(coefficients)
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
        let p = Polynomial1B7::interpolate_from_roots(&points);
        for point in points {
            assert_eq!(p.evaluate(&point), 0.into());
        }
    }

    #[test]
    fn test_interpolation_and_evaluation_zpoly() {
        let max: u32 = rand::random();
        let modulus_for_test_efficiency = 10_000;
        for point in 0..(max % modulus_for_test_efficiency) {
            assert_eq!(
                Polynomial1B7::interpolate_and_evaluate_zpoly(
                    0..(max % modulus_for_test_efficiency),
                    &PrimeFieldElement::<1_000_000_007>::from(point)
                ),
                0.into()
            );
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
        let p = Polynomial1B7::interpolate_from_coordinates(&points, &values).unwrap();
        for (x, y) in points.into_iter().zip(values) {
            assert_eq!(p.evaluate(&x), y);
        }
    }

    #[test]
    fn test_interpolation_and_evaluation_from_coordinates() {
        let number_of_points: u32 = rand::random_range(2..=100);
        let points: Vec<PrimeFieldElement<1_000_000_007>> =
            (0..number_of_points).map(PrimeFieldElement::from).collect();
        let values: Vec<PrimeFieldElement<1_000_000_007>> = (0..number_of_points)
            .map(|_| {
                let y: u32 = rand::random();
                PrimeFieldElement::from(y)
            })
            .collect();
        let p = Polynomial1B7::interpolate_from_coordinates(&points, &values).unwrap();
        for (x, y) in points.iter().zip(&values) {
            assert_eq!(
                Polynomial1B7::interpolate_and_evaluate_from_coordinates(&points, &values, x)
                    .unwrap(),
                *y
            );
        }
        let additional_points: Vec<PrimeFieldElement<1_000_000_007>> = (0..1_000)
            .map(|_| {
                let x: u32 = rand::random();
                PrimeFieldElement::from(x)
            })
            .collect();
        for x in &additional_points {
            assert_eq!(
                p.evaluate(x),
                Polynomial1B7::interpolate_and_evaluate_from_coordinates(&points, &values, x)
                    .unwrap()
            );
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
        assert_eq!(p.evaluate(&1.into()), 7.into());
        assert_eq!(p.evaluate(&0.into()), 1.into());
        let p = Polynomial1B7::new([2, 1_000_000_005, 4].map(PrimeFieldElement::from).to_vec());
        assert_eq!(p.evaluate(&0.into()), 2.into());
        assert_eq!(p.evaluate(&1.into()), 4.into());
        assert_eq!(p.evaluate(&4.into()), 58.into());
    }

    #[test]
    fn test_fft_evaluation() {
        const N: u32 = 337;
        let g = PrimeFieldElement::<N>::from(85);
        let mut units = vec![PrimeFieldElement::<N>::from(1)];
        let mut power_of_g = g;
        while power_of_g != 1.into() {
            units.push(power_of_g);
            power_of_g = power_of_g.mul(&g);
        }

        for coefficients in [
            vec![3, 1, 4, 1, 5, 9, 2, 6],
            vec![1, 0, 1, 0, 1, 0],
            vec![1, 0, 1, 0, 1],
            vec![1, 0, 0, 0],
        ] {
            let p = Polynomial::<N>::new(
                coefficients
                    .into_iter()
                    .map(PrimeFieldElement::<N>::from)
                    .collect(),
            );
            let evaluations = p.fft_evaluate(&units);
            for (i, u) in units.iter().enumerate() {
                assert_eq!(evaluations[i], p.evaluate(u));
            }
        }
    }

    #[test]
    fn test_evaluation_as_binomial() {
        // P(x) = x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        // with y = x^3 -> G(x, y) = y^2 + yx^2 + yx + y + x^2 + x + 1
        let p = Polynomial1B7::new([1, 1, 1, 1, 1, 1, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        assert_eq!(p.evaluate_as_binomial(&1.into(), &2.into(), 3), 13.into()); // 4 + 2 + 2 + 2 + 1 + 1 + 1 
        assert_eq!(p.evaluate_as_binomial(&3.into(), &2.into(), 3), 43.into()); // 4 + 18 + 6 + 2 + 9 + 3 + 1
        assert_eq!(p.evaluate_as_binomial(&2.into(), &0.into(), 3), 7.into());
        assert_eq!(p.evaluate_as_binomial(&1.into(), &1.into(), 3), 7.into());
    }

    #[test]
    fn test_partial_evaluation_as_binomial() {
        // P(x) = x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        // with y = x^3 -> G(x, y) = y^2 + yx^2 + yx + y + x^2 + x + 1
        // with x = 1 => G(1, y) = y^2 + 3y + 3
        // with x = 3 => G(3, y) = y^2 + (9 + 3 + 1)y + 9 + 3 + 1 = y^2 + 13y + 13
        let p = Polynomial1B7::new([1, 1, 1, 1, 1, 1, 1].map(PrimeFieldElement::from).to_vec()); // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
        assert_eq!(
            p.partially_evaluate_as_binomial(&1.into(), 3),
            Polynomial1B7::new(vec![3.into(), 3.into(), 1.into()])
        );
        assert_eq!(
            p.partially_evaluate_as_binomial(&3.into(), 3),
            Polynomial1B7::new(vec![13.into(), 13.into(), 1.into()])
        );

        assert_eq!(
            p.partially_evaluate_as_binomial(&3.into(), 2)
                .evaluate(&1.into()),
            p.evaluate_as_binomial(&3.into(), &1.into(), 2)
        );
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
