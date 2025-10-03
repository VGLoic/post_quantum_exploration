# Post-Quantum Exploration

This repository contains a collection of Rust implementations and experiments related to post-quantum cryptography, focusing on hash-based signature schemes and Merkle tree constructions.

This is designed for educational and experimental purposes; not optimized for production use.

## Contents

- **Winternitz One-Time Signature (WOTS):**  
  Implementation and tests for the Winternitz signature scheme, a hash-based one-time signature protocol. See [`src/winternitz_sig.rs`](src/winternitz_sig.rs).

- **Merkle Trees:**  
  Multiple versions of Merkle tree implementations for use with hash-based signatures:
  - [`src/merkletree_v0.rs`](src/merkletree_v0.rs)
  - [`src/merkletree_v1.rs`](src/merkletree_v1.rs)
  - [`src/merkletree_v2.rs`](src/merkletree_v2.rs)

- **Stateless Merkle + Winternitz Signatures:**  
  Combines Merkle trees with Winternitz signatures in a stateless fashion. See [`src/merkle_stateless_winternitz_sig.rs`](src/merkle_stateless_winternitz_sig.rs).

- **Learning with errors (LWE):**  
  A toy implementation of the Learning With Errors (LWE) problem. See [`src/toy_lwe.rs`](src/toy_lwe.rs).
  It uses vector and matrix operations over a prime field defined in [`src/primefield.rs`](src/primefield.rs).

- **Ring Learning with errors (RLWE):**  
  A toy implementation of the Ring Learning With Errors (RLWE) problem. See [`src/toy_rlwe.rs`](src/toy_rlwe.rs).
  It uses polynomial operations over a prime field defined in [`src/primefield.rs`](src/primefield.rs) and over a ring of polynomials.

- **Centered Modular Arithmetic:**  
  Implements modular arithmetic with centered coordinates for cryptographic primitives. See [`src/centered_modular_arithmetic.rs`](src/centered_modular_arithmetic.rs).

- **Prime Field Arithmetic:**  
  Provides arithmetic operations over prime fields. See [`src/primefield.rs`](src/primefield.rs).

- **NTRU Encryption:**  
  A toy implementation of the NTRU lattice-based encryption scheme. See [`src/ntru.rs`](src/ntru.rs).

- **General Modular Arithmetic:**  
  Utility functions for modular addition, multiplication, negation, and inversion. See [`src/modular_arithmetic.rs`](src/modular_arithmetic.rs).

## Getting Started

To build and test the project:

```sh
cargo build
cargo test
```

## License

This project is provided for educational purposes. See the source files for details.
