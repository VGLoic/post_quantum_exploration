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

## Getting Started

To build and test the project:

```sh
cargo build
cargo test
```

## License

This project is provided for educational purposes. See the source files for details.
