use post_quantum_exploration::{primefield::PrimeFieldElement, stark::polynomial::Polynomial};

const N: u32 = 1_000_000_007;
const P_DEGREE: usize = 10_000;
const Z_DEGREE: usize = 100_000;

type Polynomial1B7 = Polynomial<N>;
type FieldElement1B7 = PrimeFieldElement<N>;

fn main() {
    let proof = generate_proof();

    println!("Done: {proof:?}!")
}

fn generate_proof() -> Option<u32> {
    let mut points = Vec::with_capacity(P_DEGREE);
    let mut values = Vec::with_capacity(P_DEGREE);
    let mut z_roots = Vec::with_capacity(Z_DEGREE);

    for i in 0_u32..(Z_DEGREE as u32) {
        let i_as_field_element = FieldElement1B7::from(i);
        if (i as usize) < P_DEGREE {
            points.push(i_as_field_element);
            let y = rand::random_range(5u32..2_000u32);
            values.push(FieldElement1B7::from(y));
        }
        z_roots.push(i_as_field_element);
    }

    println!("Points and values have been generated");

    let _p = Polynomial1B7::interpolate_from_coordinates(points, values)?;
    println!("P polynomial has been interpolated");
    let _constraint_polynomial =
        Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
    println!("Constraint polynomial has been interpolated");
    let _z_polynomial = Polynomial1B7::interpolate_from_roots(z_roots);
    println!("Z polynomial has been interpolated");

    Some(3)
}
