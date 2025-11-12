use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use post_quantum_exploration::{
    primefield::PrimeFieldElement,
    stark::{
        polynomial::Polynomial,
        range_check::{derive_units, generate_stark_proof, verify_stark_proof},
    },
};
use std::hint::black_box;

fn bench_goldilocks(c: &mut Criterion) {
    const N: u64 = 18_446_744_069_414_584_321; // Goldilocks prime 2^64 - 2^32 + 1, see https://xn--2-umb.com/22/goldilocks/
    let mut goldilocks_group = c.benchmark_group("goldilocks");
    goldilocks_group.sample_size(20);

    for subgroup_power in [15, 18, 21].iter() {
        let generator = PrimeFieldElement::<N>::from(7);
        let subgroup_exponent = (N - 1) / (1 << subgroup_power);
        let subgroup_generator = generator.exp(subgroup_exponent);
        goldilocks_group.bench_with_input(
            BenchmarkId::new("units derivation", subgroup_power),
            subgroup_power,
            |b, _| {
                b.iter(|| {
                    derive_units(black_box(subgroup_generator));
                })
            },
        );

        let p_max_degree = 1u64 << (subgroup_power - 10);

        let units = derive_units(subgroup_generator);
        let mut polynomial = Polynomial::<N>::interpolate_from_roots(
            (1..p_max_degree)
                .map(PrimeFieldElement::<N>::from)
                .collect::<Vec<PrimeFieldElement<N>>>()
                .as_slice(),
        );
        polynomial =
            polynomial.mul_by_scalar(&polynomial.evaluate(&p_max_degree.into()).inv().unwrap());

        goldilocks_group.bench_with_input(
            BenchmarkId::new("fft", subgroup_power),
            subgroup_power,
            |b, _| {
                b.iter(|| {
                    polynomial.fft_evaluate(black_box(&units));
                });
            },
        );

        goldilocks_group.bench_with_input(
            BenchmarkId::new("stark proof generation", subgroup_power),
            subgroup_power,
            |b, _| {
                b.iter(|| {
                    generate_stark_proof(
                        black_box(polynomial.clone()),
                        black_box(subgroup_generator),
                        black_box(p_max_degree),
                    )
                    .unwrap();
                });
            },
        );

        let proof = generate_stark_proof(polynomial, subgroup_generator, p_max_degree).unwrap();

        goldilocks_group.bench_with_input(
            BenchmarkId::new("stark proof verification", subgroup_power),
            subgroup_power,
            |b, _| {
                b.iter(|| {
                    verify_stark_proof(
                        black_box(proof.clone()),
                        black_box(subgroup_generator),
                        black_box(p_max_degree),
                    )
                    .unwrap();
                });
            },
        );
    }

    goldilocks_group.finish();
}

criterion_group!(benches, bench_goldilocks);
criterion_main!(benches);
