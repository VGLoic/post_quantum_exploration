#[cfg(test)]
mod merkletree_tests {
    use bitvec::prelude::*;
    use std::num::TryFromIntError;

    use anyhow::anyhow;
    use sha3::{Digest, Sha3_256};

    // ##########################################################
    // ##########################################################
    // ##########################################################
    // ###################### V1 ################################
    // ##########################################################
    // ##########################################################
    // ##########################################################

    /*
     * # Construction of Merkle tree v1
     * Inputs:
     * - list of values,
     * - depth
     * Algorithm:
     * 1. resize the list of values according to depth,
     * 2. create the leaves,
     * 3. creates the first level of branches with pair of leaves,
     * 4. creates new levels of branches with pair of previously constructed branches,
     * 5. there should be a single branch once iterated d times,
     * 6. build the tree with it.
     *
     * Complexity wise:
     * - O(N) for creating the values,
     * - then we do `log(N)` loops of order O(log(N)),
     * - the complexity is therefore O(N)
     */

    #[derive(Clone)]
    struct LeafV1 {
        value: [u8; 32],
        hash: [u8; 32],
    }
    struct BranchV1 {
        left: Box<NodeV1>,
        right: Box<NodeV1>,
        hash: [u8; 32],
    }

    enum NodeV1 {
        Leaf(LeafV1),
        Branch(BranchV1),
    }

    impl NodeV1 {
        fn hash(&self) -> &[u8; 32] {
            match self {
                Self::Branch(b) => &b.hash,
                Self::Leaf(l) => &l.hash,
            }
        }
    }

    struct MerkleTreeV1 {
        depth: u8,
        root: BranchV1,
    }

    struct ValueWithProof<'a> {
        value: &'a [u8; 32],
        proof: Vec<(bool, &'a [u8; 32])>,
    }

    impl MerkleTreeV1 {
        fn root_hash(&self) -> &[u8; 32] {
            &self.root.hash
        }

        fn get<'a>(&'a self, selector: &[bool]) -> Result<ValueWithProof<'a>, anyhow::Error> {
            if selector.len() != self.depth.into() {
                return Err(anyhow!(
                    "invalid selector, expected of length: {}, got length {}",
                    self.depth,
                    selector.len()
                ));
            }
            let mut current_branch: &BranchV1 = &self.root;
            let mut proof: Vec<(bool, &'a [u8; 32])> = vec![];
            for direction in selector {
                // If direction is true, we register the hash of the left node in the proof and we take a look at the right node

                let next_node = if *direction {
                    proof.push((false, current_branch.left.hash()));
                    &current_branch.right
                } else {
                    // If direction is false, we register the hash of the right node in the proof and we take a look at the left node
                    proof.push((true, current_branch.right.hash()));
                    &current_branch.left
                };
                match &**next_node {
                    NodeV1::Branch(b) => {
                        current_branch = b;
                    }
                    NodeV1::Leaf(l) => {
                        proof.reverse();
                        return Ok(ValueWithProof {
                            value: &l.value,
                            proof,
                        });
                    }
                };
            }

            Err(anyhow!(
                "unexpected end of direction selector without leaves"
            ))
        }

        fn new(depth: u8, values: &[[u8; 32]]) -> Result<MerkleTreeV1, anyhow::Error> {
            if depth == 0 {
                return Err(anyhow!("depth must be more than 0"));
            }
            if depth > 32 {
                return Err(anyhow!("depth must be less than 32, got {depth}"));
            }
            let number_of_leaves = (2u32).pow(depth.into());
            let number_of_leaves_usize =
                number_of_leaves.try_into().map_err(|e: TryFromIntError| {
                    anyhow!(e).context("failed to convert number of leaves to usize")
                })?;
            let number_input_values = values.len();
            if number_input_values > number_of_leaves_usize {
                return Err(anyhow!("too many values, do not fill in tree"));
            }

            // We build all the leaves
            let mut leaves: Vec<LeafV1> = values
                .iter()
                .map(|v| LeafV1 {
                    value: v.to_owned(),
                    hash: Sha3_256::digest(v).into(),
                })
                .collect();

            if number_input_values < number_of_leaves_usize {
                let zero_value = [0u8; 32];
                let zero_hash: [u8; 32] = Sha3_256::digest(zero_value).into();
                leaves.resize(
                    number_of_leaves_usize,
                    LeafV1 {
                        value: zero_value,
                        hash: zero_hash,
                    },
                );
            }

            // We pair leaves in order to build the branches of level `d - 1`
            let mut branches: Vec<Box<NodeV1>> = vec![];
            let mut pair: Vec<LeafV1> = vec![];
            for leaf in leaves.into_iter() {
                pair.push(leaf);
                if pair.len() == 2 {
                    let right = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                    let left = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                    let mut hasher = Sha3_256::new();
                    hasher.update(left.hash);
                    hasher.update(right.hash);
                    let branch = BranchV1 {
                        left: Box::new(NodeV1::Leaf(left)),
                        right: Box::new(NodeV1::Leaf(right)),
                        hash: hasher.finalize().into(),
                    };
                    let node = NodeV1::Branch(branch);
                    branches.push(Box::new(node));
                    pair.clear();
                }
            }

            // We iteratively pair branches/nodes until we reach depth
            for _ in 1..depth {
                let mut next_branches: Vec<Box<NodeV1>> = vec![];
                let mut pair: Vec<Box<NodeV1>> = vec![];
                for b in branches.into_iter() {
                    pair.push(b);
                    if pair.len() == 2 {
                        let right = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                        let left = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                        let mut hasher = Sha3_256::new();
                        hasher.update(left.hash());
                        hasher.update(right.hash());
                        let branch = BranchV1 {
                            left,
                            right,
                            hash: hasher.finalize().into(),
                        };
                        let node = NodeV1::Branch(branch);
                        next_branches.push(Box::new(node));
                        pair.clear();
                    }
                }
                branches = next_branches;
            }

            // At the end of iteration, we should be left with exactly one branch
            if branches.len() != 1 {
                return Err(anyhow!(
                    "found more branches than expected! Got {} instead of 1",
                    branches.len()
                ));
            }
            let root = match *branches
                .pop()
                .ok_or(anyhow!("branches should have at least one element"))?
            {
                NodeV1::Branch(b) => b,
                NodeV1::Leaf(_) => return Err(anyhow!("unexpected leaf at the root retrieval")),
            };

            Ok(MerkleTreeV1 { depth, root })
        }
    }

    fn verify_proof(
        root: &[u8; 32],
        value: &[u8; 32],
        proof: &[(bool, &[u8; 32])],
    ) -> Result<(), anyhow::Error> {
        let mut recovered_hash = Sha3_256::digest(value);
        // We iterate over the proof elements, if first element is true, it means the accumulator should be hashed first, else the accumulator should be hashed second
        for (direction, sibling_hash) in proof {
            let mut hasher = Sha3_256::new();
            if *direction {
                hasher.update(recovered_hash);
                hasher.update(sibling_hash);
            } else {
                hasher.update(sibling_hash);
                hasher.update(recovered_hash);
            }
            recovered_hash = hasher.finalize();
        }
        if &Into::<[u8; 32]>::into(recovered_hash) != root {
            return Err(anyhow!("invalid proof"));
        }

        Ok(())
    }

    #[test]
    fn test_v1() {
        let depth = 16u8;
        let number_of_elements = 2usize.pow(depth.into());
        let mut values = vec![];
        for _ in 0..number_of_elements {
            let inner_value: u8 = rand::random();
            values.push([inner_value; 32]);
        }
        let tree = MerkleTreeV1::new(depth, &values).unwrap();

        for (i, value) in values.iter().enumerate() {
            let selector = i
                .view_bits::<Lsb0>()
                .into_iter()
                .take(depth.into())
                .map(|v| *v)
                .rev()
                .collect::<Vec<bool>>();
            let selected = tree.get(&selector).unwrap();
            assert_eq!(selected.value, value);
            if let Err(e) =
                verify_proof(tree.root_hash(), selected.value, selected.proof.as_slice())
            {
                panic!("{e:?}");
            }
        }
    }
}
