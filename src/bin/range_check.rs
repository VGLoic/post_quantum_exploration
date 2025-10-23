use post_quantum_exploration::{
    primefield::PrimeFieldElement,
    stark::{
        polynomial::Polynomial,
        range_check::{generate_stark_proof, verify_stark_proof},
    },
};

const N: u32 = 1_073_153; // Chosen because: (p - 1) is divisible 6 times by 4 and is larger than 1_000_000
const GENERATOR: u32 = 3;
const P_MAX_DEGREE: u32 = 1_024; // 0 <= P(x) <= 9 for 1 <= x <= P_MAX_DEGREE

fn main() {
    let mut p = Polynomial::<N>::interpolate_from_roots(
        (1..P_MAX_DEGREE)
            .map(PrimeFieldElement::<N>::from)
            .collect::<Vec<PrimeFieldElement<N>>>()
            .as_slice(),
    );
    p = p.mul_by_scalar(&p.evaluate(&P_MAX_DEGREE.into()).inv().unwrap());

    println!("Start generating stark proof");

    let stark_proof = generate_stark_proof(p, GENERATOR.into(), P_MAX_DEGREE).unwrap();

    println!("Stark proof successfully generated");

    verify_stark_proof(stark_proof, GENERATOR.into(), P_MAX_DEGREE).unwrap();

    println!("Stark proof successfully verified");
}
