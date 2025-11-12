use anyhow::anyhow;

use crate::primefield::PrimeFieldElement;

#[allow(dead_code)]
/// Performs the Fourier Transform over a whole set of units generated as successive power of a generator.
///
/// The implementation uses the fast Fourier transform, therefore the units set must have a length as 2^k.
///
/// The number of values must be equal to the length of the units set.
///
/// # Arguments
/// * `values` - Set of values on which the Fourier transform is applied. The set must have the same length as the units,
/// * `units` - Set of units, generated as successive power of a generator. The set must have a length of 2^k.
///
/// Returns the Fourier transform of the values over the whole set.
pub fn fft_recursive<const N: u64>(
    values: &[PrimeFieldElement<N>],
    units: &[PrimeFieldElement<N>],
) -> Vec<PrimeFieldElement<N>> {
    if values.len() == 1 {
        return vec![values[0]; units.len()];
    }

    let mut odd_values = Vec::with_capacity(values.len() / 2);
    let mut even_values = Vec::with_capacity(values.len() / 2);

    for (i, &v) in values.iter().enumerate() {
        if i & 1 == 0 {
            even_values.push(v);
        } else {
            odd_values.push(v);
        }
    }

    let reduced_units: Vec<PrimeFieldElement<N>> = units.iter().step_by(2).copied().collect();
    let even_values = fft(&even_values, &reduced_units);
    let odd_values = fft(&odd_values, &reduced_units);

    let zero = PrimeFieldElement::<N>::from(0);
    let mut result = vec![zero; units.len()];
    let reduced_units_len = units.len() / 2;
    for i in 0..reduced_units_len {
        let odd_contribution = units[i].mul(&odd_values[i]);
        result[i] = even_values[i].add(&odd_contribution);
        result[i + reduced_units_len] = even_values[i].add(&odd_contribution.neg());
    }
    result
}

pub fn fft<const N: u64>(
    values: &[PrimeFieldElement<N>],
    units: &[PrimeFieldElement<N>],
) -> Vec<PrimeFieldElement<N>> {
    let reductions = values.len().trailing_zeros();
    let mut result = reorganize_fourier_values(values, reductions);

    for reduction_level in 0..reductions {
        let batch_size = 1 << reduction_level; // 1, 2, 4, 8, ..
        let units_step = 1 << (reductions - reduction_level - 1); // ..., 8, 4, 2, 1,
        // At iteration 0, we take 1 element on 2, e.g. 0, 2, 4, 6, 8, 10, 12, 14, ...
        // At iteration 1, we take 2 elements, skip 2, etc... e.g. 0, 1, 4, 5, 8, 9, 12, 13, ...
        // At iteration 2, we take 0, 1, 2, 3, 8, 9, 10, 11, ...
        let mut iteration_count = 0;
        let mut i = 0;
        while i < units.len() {
            if iteration_count == batch_size {
                iteration_count = 0;
                i += batch_size;
            } else {
                let conjugated_index = i + batch_size;
                let unit = &units[units_step * (i % batch_size)];
                let odd_contribution = unit.mul(&result[conjugated_index]);
                (result[i], result[conjugated_index]) = (
                    result[i].add(&odd_contribution),
                    result[i].add(&odd_contribution.neg()),
                );
                iteration_count += 1;
                i += 1;
            }
        }
    }

    result
}

fn reorganize_fourier_values<const N: u64>(
    values: &[PrimeFieldElement<N>],
    reductions: u32,
) -> Vec<PrimeFieldElement<N>> {
    let mut result = Vec::with_capacity(values.len());
    let mut indices = Vec::with_capacity(values.len());
    indices.push(0);
    result.push(values[0]);
    for reduction_level in 1..=reductions {
        let lead_power_of_two = 1 << (reductions - reduction_level);
        for j in 0..(1 << (reduction_level - 1)) {
            let index = lead_power_of_two + indices[j];
            indices.push(index);
            result.push(values[index]);
        }
    }
    result
}

#[allow(dead_code)]
pub fn inverse_fft<const N: u64>(
    values: &[PrimeFieldElement<N>],
    units: &[PrimeFieldElement<N>],
) -> Result<Vec<PrimeFieldElement<N>>, anyhow::Error> {
    if values.is_empty() {
        return Ok(vec![0.into(); units.len()]);
    }
    let ft = fft(values, units);
    let len_as_u64: u64 = units
        .len()
        .try_into()
        .map_err(|_| anyhow!("units are too long, can not express the length as u32"))?;
    let inv_len = PrimeFieldElement::<N>::from(len_as_u64)
        .inv()
        .ok_or(anyhow!("unable to inverse the length of units"))?;

    let mut result = Vec::with_capacity(units.len());

    let first_value = ft.first().unwrap_or_else(|| {
        unreachable!("there has been a check on the emptyness of values before")
    });
    result.push(first_value.mul(&inv_len));

    for fv in ft.iter().skip(1).rev() {
        result.push(fv.mul(&inv_len));
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::{
        primefield::PrimeFieldElement,
        stark::fft::{fft, fft_recursive, inverse_fft},
    };

    #[test]
    fn test_fft() {
        const N: u64 = 337;
        let g = PrimeFieldElement::<N>::from(85);
        let mut units = vec![PrimeFieldElement::<N>::from(1)];
        let mut power_of_g = g;
        while power_of_g != 1.into() {
            units.push(power_of_g);
            power_of_g = power_of_g.mul(&g);
        }

        for values in [
            vec![3, 1, 4, 1, 5, 9, 2, 6],
            vec![1, 0, 1, 0, 1, 0],
            vec![1, 0, 1, 0, 1],
            vec![1u64, 0, 0, 0],
        ] {
            let mut values: Vec<PrimeFieldElement<N>> = values
                .iter()
                .map(|v| PrimeFieldElement::<N>::from(*v))
                .collect();
            values.resize(units.len(), 0.into());
            let ft = fft(&values, &units);
            let ft_recursive = fft_recursive(&values, &units);
            assert_eq!(ft, ft_recursive);

            assert_eq!(inverse_fft(&ft, &units).unwrap(), values);
        }
    }
}
