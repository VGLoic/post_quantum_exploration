/*
 * Toy example of stateless Winternitz signatures using
 * - message size of 16 bits,
 * - Winternitz parameters: w = 16 = 2^4
 * - we use a Merkle tree of depth 16 to handle the 2^16 possible keys,
 *
 *
 * We end up with groups of 4 bits, hence y = 16 / 4 = 4 and checksum is at most 4 * 16 = 64 = 2^6
 * Therefore: 6 (4 + 2) signatures and public keys
 */
