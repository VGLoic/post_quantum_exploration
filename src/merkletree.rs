#[cfg(test)]
mod merkletree_tests {
    use bitvec::prelude::*;
    use std::{num::TryFromIntError, slice::Iter};

    use anyhow::anyhow;
    use sha3::{Digest, Sha3_256};

    // ##########################################################
    // ##########################################################
    // ##########################################################
    // ###################### V0 ################################
    // ##########################################################
    // ##########################################################
    // ##########################################################

    struct LeafV0 {
        value: [u8; 32],
        hash: [u8; 32],
    }
    struct BranchV0 {
        left_node: Option<Box<NodeV0>>,
        right_node: Option<Box<NodeV0>>,
        hash: Option<[u8; 32]>,
    }

    enum NodeV0 {
        Leaf(LeafV0),
        Branch(BranchV0),
    }

    struct MerkleTreeV0 {
        depth: u8,
        root: BranchV0,
    }

    // ##########################################################
    // ###################### CONSTRUCTION ######################
    // ##########################################################

    /*
     * # Construction of Merkle tree v0
     * Inputs:
     * - list of values,
     * - depth
     * Algorithm:
     * 1. resize the list of values according to depth,
     * 2. create the root with empty left & right
     * 3. iterate over the values, starting from the root:
     *  3.1. if left is empty,
     *      3.1.a. if we are at depth, creates a leaf and SUCCESS
     *      3.1.b. else, creates an empty branch, link it to parent and restart the loop,
     *  3.2. if left is not empty,
     *      3.2.a. if node is a leaf, WE CAN'T ADD IT HERE, WE ABORT,
     *      3.2.b. else, node is branch, we restart the loop from it,
     *  3.3. if right is empty,
     *      3.3.a. if we are at depth, creates a leaf and SUCCESS
     *      3.3.b. else, creates an empty branch, link it to parent and restart the loop,
     *  3.4. if right is not empty,
     *      3.4.a. if node is a leaf, WE CAN'T ADD IT HERE, WE ABORT,
     *      3.4.b else, node is a branch, we restart the loop from it.
     *
     *
     * Complexity:
     *  - iteration over the values: 2^d = O(N)
     *  - for each value, we check the whole tree. For last value we check 2^d + 2^(d-1) + 2^(d-2) + ... = O(N)
     *  - so we get O(N^2)
     */

    impl NodeV0 {
        fn try_hash(&self) -> Option<&[u8; 32]> {
            match self {
                Self::Branch(b) => b.hash.as_ref(),
                Self::Leaf(l) => Some(&l.hash),
            }
        }
    }

    impl BranchV0 {
        fn try_to_add(
            &mut self,
            v: [u8; 32],
            current_depth: u8,
            depth: u8,
        ) -> Result<(), anyhow::Error> {
            if self.try_to_add_left(v, current_depth, depth).is_err() {
                self.try_to_add_right(v, current_depth, depth)?;
                if let Some(left) = &self.left_node
                    && let Some(right) = &self.right_node
                    && let Some(left_hash) = left.try_hash()
                    && let Some(right_hash) = right.try_hash()
                {
                    let mut hasher = Sha3_256::new();
                    hasher.update(left_hash);
                    hasher.update(right_hash);
                    let hash = hasher.finalize();
                    self.hash = Some(hash.into());
                }
            }
            Ok(())
        }

        fn try_to_add_left(
            &mut self,
            v: [u8; 32],
            current_depth: u8,
            depth: u8,
        ) -> Result<(), anyhow::Error> {
            match &mut self.left_node {
                None => {
                    // If we are at the bottom, we add a new leaf
                    if current_depth == depth - 1 {
                        let leaf = LeafV0 {
                            value: v,
                            hash: Sha3_256::digest(v).into(),
                        };
                        self.left_node = Some(Box::new(NodeV0::Leaf(leaf)));
                        Ok(())
                    } else {
                        // Else, we build a new branch and we try to add the value to it by increasing the current depth
                        let mut branch = BranchV0 {
                            left_node: None,
                            right_node: None,
                            hash: None,
                        };
                        branch.try_to_add(v, current_depth + 1, depth)?;
                        self.left_node = Some(Box::new(NodeV0::Branch(branch)));
                        Ok(())
                    }
                }
                Some(node) => match &mut **node {
                    NodeV0::Leaf(_) => Err(anyhow!("can not add here")),
                    NodeV0::Branch(b) => b.try_to_add(v, current_depth + 1, depth),
                },
            }
        }

        fn try_to_add_right(
            &mut self,
            v: [u8; 32],
            current_depth: u8,
            depth: u8,
        ) -> Result<(), anyhow::Error> {
            match &mut self.right_node {
                None => {
                    // If we are at the bottom, we add a new leaf
                    if current_depth == depth - 1 {
                        let leaf = LeafV0 {
                            value: v,
                            hash: Sha3_256::digest(v).into(),
                        };
                        self.right_node = Some(Box::new(NodeV0::Leaf(leaf)));
                        Ok(())
                    } else {
                        // Else, we build a new branch and we try to add the value to it by increasing the current depth
                        let mut branch = BranchV0 {
                            left_node: None,
                            right_node: None,
                            hash: None,
                        };
                        branch.try_to_add(v, current_depth + 1, depth)?;
                        self.right_node = Some(Box::new(NodeV0::Branch(branch)));
                        Ok(())
                    }
                }
                Some(node) => match &mut **node {
                    NodeV0::Leaf(_) => Err(anyhow!("can not add here")),
                    NodeV0::Branch(b) => b.try_to_add(v, current_depth + 1, depth),
                },
            }
        }
    }

    impl MerkleTreeV0 {
        fn add(&mut self, v: [u8; 32]) -> Result<(), anyhow::Error> {
            if self.root.try_to_add_left(v, 0, self.depth).is_err()
                && self.root.try_to_add_right(v, 0, self.depth).is_err()
            {
                return Err(anyhow!("tree is full"));
            }

            Ok(())
        }

        fn new(depth: u8, values: &[[u8; 32]]) -> Result<MerkleTreeV0, anyhow::Error> {
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

            let root = BranchV0 {
                left_node: None,
                right_node: None,
                hash: None,
            };

            let mut tree = MerkleTreeV0 { root, depth };

            for v in values {
                tree.add(v.to_owned())?;
            }

            let remaining_values_to_add = number_of_leaves_usize - number_input_values;

            for _ in 0..remaining_values_to_add {
                tree.add([0u8; 32])?;
            }

            if let Some(left) = &tree.root.left_node
                && let Some(left_hash) = left.try_hash()
                && let Some(right) = &tree.root.right_node
                && let Some(right_hash) = right.try_hash()
            {
                let mut hasher = Sha3_256::new();
                hasher.update(left_hash);
                hasher.update(right_hash);
                tree.root.hash = Some(hasher.finalize().into());
            } else {
                return Err(anyhow!("filling the tree was wrong"));
            }

            Ok(tree)
        }
    }

    // #######################################################
    // ###################### RETRIEVAL ######################
    // #######################################################

    struct ValueWithProof<'a> {
        value: &'a [u8; 32],
        proof: Vec<(bool, &'a [u8; 32])>,
    }

    impl NodeV0 {
        fn hash(&self) -> &[u8; 32] {
            match self {
                Self::Branch(b) => b.hash.as_ref().unwrap_or_else(|| {
                    unreachable!("hash of a branch is always defined once initialized")
                }),
                Self::Leaf(l) => &l.hash,
            }
        }

        fn get<'a>(&'a self, mut selector: Iter<'a, bool>) -> ValueWithProof<'a> {
            match self {
                Self::Leaf(l) => ValueWithProof {
                    value: &l.value,
                    proof: vec![],
                },
                Self::Branch(b) => {
                    let s = selector.next().unwrap_or_else(|| {
                        unreachable!("iterator always contains the proper number of elements")
                    });
                    if *s {
                        let child_value = b.right().get(selector);
                        let mut proof = child_value.proof;
                        proof.push((false, b.left().hash()));
                        ValueWithProof {
                            value: child_value.value,
                            proof,
                        }
                    } else {
                        let child_value = b.left().get(selector);
                        let mut proof = child_value.proof;
                        proof.push((true, b.right().hash()));
                        ValueWithProof {
                            value: child_value.value,
                            proof,
                        }
                    }
                }
            }
        }
    }

    impl BranchV0 {
        fn left(&self) -> &NodeV0 {
            self.left_node
                .as_ref()
                .unwrap_or_else(|| unreachable!("left of a branch is always defined"))
        }
        fn right(&self) -> &NodeV0 {
            self.right_node
                .as_ref()
                .unwrap_or_else(|| unreachable!("right of a branch is always defined"))
        }
    }

    impl MerkleTreeV0 {
        fn left(&self) -> &NodeV0 {
            self.root.left()
        }

        fn right(&self) -> &NodeV0 {
            self.root.right()
        }

        fn get<'a>(&'a self, selector: &'a [bool]) -> Result<ValueWithProof<'a>, anyhow::Error> {
            if selector.len() != self.depth.into() {
                return Err(anyhow!("selector not long enough"));
            }

            let mut iterator = selector.iter();
            let first = iterator
                .next()
                .unwrap_or_else(|| unreachable!("iterator can not be empty"));
            if *first {
                let child_value = self.right().get(iterator);
                let mut proof = child_value.proof;
                proof.push((false, self.left().hash()));
                Ok(ValueWithProof {
                    value: child_value.value,
                    proof,
                })
            } else {
                let child_value = self.left().get(iterator);
                let mut proof = child_value.proof;
                proof.push((true, self.right().hash()));
                Ok(ValueWithProof {
                    value: child_value.value,
                    proof,
                })
            }
        }

        fn verify_proof(
            &self,
            value: &[u8; 32],
            proof: &[(bool, &[u8; 32])],
        ) -> Result<(), anyhow::Error> {
            let mut recovered_root_hash: [u8; 32] = Sha3_256::digest(value).into();
            for (direction, sibling_hash) in proof {
                let mut hasher = Sha3_256::new();
                if *direction {
                    hasher.update(recovered_root_hash);
                    hasher.update(sibling_hash);
                } else {
                    hasher.update(sibling_hash);
                    hasher.update(recovered_root_hash);
                }
                recovered_root_hash = hasher.finalize().into();
            }
            if recovered_root_hash
                != self
                    .root
                    .hash
                    .ok_or_else(|| anyhow!("tree not initialized"))?
            {
                return Err(anyhow!("proof is invalid"));
            }

            Ok(())
        }
    }

    #[test]
    fn test_v0() {
        let depth = 13u8;
        let number_of_elements = 2usize.pow(depth.into());
        let mut values = vec![];
        for _ in 0..number_of_elements {
            let inner_value: u8 = rand::random();
            values.push([inner_value; 32]);
        }
        let tree = MerkleTreeV0::new(depth, &values).unwrap();

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
            if let Err(e) = tree.verify_proof(selected.value, selected.proof.as_slice()) {
                panic!("{e:?}");
            }
        }
    }

    // ##########################################################
    // ##########################################################
    // ##########################################################
    // ###################### V1 ################################
    // ##########################################################
    // ##########################################################
    // ##########################################################

    // ##########################################################
    // ###################### CONSTRUCTION ######################
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

    impl MerkleTreeV1 {
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

        fn verify_proof(
            &self,
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
            if Into::<[u8; 32]>::into(recovered_hash) != self.root.hash {
                return Err(anyhow!("invalid proof"));
            }

            Ok(())
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

    #[test]
    fn test_v1() {
        let depth = 13u8;
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
            if let Err(e) = tree.verify_proof(selected.value, selected.proof.as_slice()) {
                panic!("{e:?}");
            }
        }
    }
}
