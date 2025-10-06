use post_quantum_exploration::stark::range_check::generate_proof;

fn main() {
    let proof = generate_proof();

    println!("Done: {proof:?}!")
}
