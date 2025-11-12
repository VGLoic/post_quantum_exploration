use post_quantum_exploration::{
    primefield::PrimeFieldElement,
    stark::{
        fri::derive_fri_reductions_count,
        polynomial::Polynomial,
        range_check::{generate_stark_proof, verify_stark_proof},
    },
};

const N: u64 = 18_446_744_069_414_584_321; // Goldilocks prime 2^64 - 2^32 + 1, see https://xn--2-umb.com/22/goldilocks/

fn main() {
    let subgroup_power: u32 = std::env::var("SUBGROUP_POWER")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(28);

    let p_max_degree: u64 = std::env::var("P_MAX_DEGREE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1 << (subgroup_power - 10));

    println!(
        "Using Goldilocks field with subgroup power {} and polynomial max degree {}",
        subgroup_power, p_max_degree
    );

    if derive_fri_reductions_count(p_max_degree).is_err() {
        panic!("P_MAX_DEGREE must be divisible by 4 enough times");
    }

    let mut p = Polynomial::<N>::interpolate_from_roots(
        (1..p_max_degree)
            .map(PrimeFieldElement::<N>::from)
            .collect::<Vec<PrimeFieldElement<N>>>()
            .as_slice(),
    );
    p = p.mul_by_scalar(&p.evaluate(&p_max_degree.into()).inv().unwrap());

    let generator = PrimeFieldElement::<N>::from(7);
    let subgroup_exponent = (N - 1) / (1 << subgroup_power);
    let subgroup_generator = generator.exp(subgroup_exponent);

    println!("Start generating stark proof");

    let stark_proof = generate_stark_proof(p, subgroup_generator, p_max_degree).unwrap();

    println!("Stark proof successfully generated");

    verify_stark_proof(stark_proof, subgroup_generator, p_max_degree).unwrap();

    println!("Stark proof successfully verified");
}
