use anyhow::anyhow;
use post_quantum_exploration::{
    merkletree_v2::MerkleTreeV2,
    primefield::PrimeFieldElement,
    stark::{polynomial::Polynomial, range_check::StarkLeaf},
};

const N: u32 = 1_000_000_007;
const P_DEGREE: usize = 100;
const Z_DEGREE: usize = 1_000;
const TOTAL_POINTS: u32 = 1_000_000;

type Polynomial1B7 = Polynomial<N>;
type FieldElement1B7 = PrimeFieldElement<N>;

fn main() {
    let proof = generate_proof();

    println!("Done: {proof:?}!")
}

fn generate_proof() -> Result<u32, anyhow::Error> {
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

    let p = Polynomial1B7::interpolate_from_coordinates(points, values).ok_or(anyhow!(
        "Unable to interpolate P polynomial from coordinates"
    ))?;
    println!("P polynomial has been interpolated");
    let constraint_polynomial =
        Polynomial1B7::interpolate_from_roots((0..=9).map(FieldElement1B7::from).collect());
    println!("Constraint polynomial has been interpolated");
    let z_polynomial = Polynomial1B7::interpolate_from_roots(z_roots);
    println!("Z polynomial has been interpolated");

    let mut values: Vec<StarkLeaf<N>> = Vec::with_capacity(1 + Z_DEGREE);
    for i in 0u32..TOTAL_POINTS {
        if i > 0 && i % 100_000 == 0 {
            println!("evaluation reached a 100_000 step");
        }
        let i_as_field_element = FieldElement1B7::from(i);

        let p_eval = p.evaluate(i_as_field_element);
        let cp_eval = constraint_polynomial.evaluate(p_eval);
        let d_eval = if (i as usize) < Z_DEGREE {
            0.into()
        } else {
            let z_eval = z_polynomial.evaluate(i_as_field_element);
            cp_eval.mul(
                &z_eval
                    .inv()
                    .ok_or(anyhow!("unable to find inverse of Z evaluation"))?,
            )
        };

        values.push(StarkLeaf::new(cp_eval, d_eval));
    }

    let _tree = MerkleTreeV2::new(20, &values)?;

    Ok(3)
}
