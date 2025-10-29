use crate::primefield::PrimeFieldElement;

/// Performs the Fourier Transform over a whole set of units generated as successive power of a generator.
///
/// The implementation uses the fast Fourier transform, therefore the units set must have a length as 2^k.
///
/// The number of values must be less or equal to the length of the units set.
///
/// # Arguments
/// * `values` - Set of values on which the Fourier transform is applied
/// * `units` - Set of units, generated as successive power of a generator. The set must have a length of 2^k.
///
/// Returns the Fourier transform of the values over the whole set.
pub fn fft<const N: u32>(
    values: &[PrimeFieldElement<N>],
    units: &[PrimeFieldElement<N>],
) -> Vec<PrimeFieldElement<N>> {
    if values.is_empty() {
        return vec![0.into(); units.len()];
    }
    if values.len() == 1 {
        return vec![values[0]; units.len()];
    }

    let mut odd_values = vec![];
    let mut even_values = vec![];

    for (i, &v) in values.iter().enumerate() {
        if i.is_multiple_of(2) {
            even_values.push(v);
        } else {
            odd_values.push(v);
        }
    }
    let zero = PrimeFieldElement::<N>::from(0);
    while let Some(v) = even_values.last()
        && v == &zero
    {
        even_values.pop();
    }
    while let Some(v) = odd_values.last()
        && v == &zero
    {
        odd_values.pop();
    }

    let reduced_units: Vec<PrimeFieldElement<N>> =
        units.iter().step_by(2).map(|v| v.to_owned()).collect();
    let even_values = fft(&even_values, &reduced_units);
    let odd_values = fft(&odd_values, &reduced_units);

    let mut first_half_result = Vec::with_capacity(units.len() / 2);
    let mut second_half_result = Vec::with_capacity(units.len() / 2);

    for i in 0..(units.len() / 2) {
        let odd_contribution = units[i].mul(&odd_values[i]);
        first_half_result.push(even_values[i].add(&odd_contribution));
        second_half_result.push(even_values[i].add(&odd_contribution.neg()));
    }

    [first_half_result, second_half_result].concat()
}
