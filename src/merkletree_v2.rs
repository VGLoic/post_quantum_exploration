use std::num::TryFromIntError;

use anyhow::anyhow;
use sha3::{Digest, Sha3_256};

// ##########################################################
// ##########################################################
// ##########################################################
// ###################### V2 ################################
// ##########################################################
// ##########################################################
// ##########################################################

/*
 * Same than V1 but with arbitrary leaf value
 */

#[allow(dead_code)]
#[derive(Clone)]
struct LeafV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    value: T,
    hash: [u8; 32],
}

#[allow(dead_code)]
struct BranchV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    left: Box<NodeV2<T>>,
    right: Box<NodeV2<T>>,
    hash: [u8; 32],
}

#[allow(dead_code)]
enum NodeV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    Leaf(LeafV2<T>),
    Branch(BranchV2<T>),
}

impl<T> NodeV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    #[allow(dead_code)]
    fn hash(&self) -> &[u8; 32] {
        match self {
            Self::Branch(b) => &b.hash,
            Self::Leaf(l) => &l.hash,
        }
    }
}

#[allow(dead_code)]
pub struct MerkleTreeV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    depth: u8,
    root: BranchV2<T>,
}

#[allow(dead_code)]
pub struct ValueWithProof<'a, T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    pub value: &'a T,
    pub proof: Vec<(bool, &'a [u8; 32])>,
}

impl<T> MerkleTreeV2<T>
where
    T: Clone + Default + AsRef<[u8]>,
{
    #[allow(dead_code)]
    pub fn root_hash(&self) -> &[u8; 32] {
        &self.root.hash
    }

    #[allow(dead_code)]
    pub fn get<'a>(&'a self, selector: &[bool]) -> Result<ValueWithProof<'a, T>, anyhow::Error> {
        if selector.len() != self.depth.into() {
            return Err(anyhow!(
                "invalid selector, expected of length: {}, got length {}",
                self.depth,
                selector.len()
            ));
        }
        let mut current_branch: &BranchV2<T> = &self.root;
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
                NodeV2::Branch(b) => {
                    current_branch = b;
                }
                NodeV2::Leaf(l) => {
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

    #[allow(dead_code)]
    pub fn new(depth: u8, values: &[T]) -> Result<MerkleTreeV2<T>, anyhow::Error> {
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
        let mut leaves: Vec<LeafV2<T>> = values
            .iter()
            .map(|v| LeafV2 {
                value: v.to_owned(),
                hash: Sha3_256::digest(v).into(),
            })
            .collect();

        if number_input_values < number_of_leaves_usize {
            let zero_value = T::default();
            let zero_hash = Sha3_256::digest(&zero_value);
            leaves.resize(
                number_of_leaves_usize,
                LeafV2 {
                    value: zero_value,
                    hash: zero_hash.into(),
                },
            );
        }

        // We pair leaves in order to build the branches of level `d - 1`
        let mut branches: Vec<Box<NodeV2<T>>> = vec![];
        let mut pair: Vec<LeafV2<T>> = vec![];
        for leaf in leaves.into_iter() {
            pair.push(leaf);
            if pair.len() == 2 {
                let right = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                let left = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                let mut hasher = Sha3_256::new();
                hasher.update(left.hash);
                hasher.update(right.hash);
                let branch = BranchV2 {
                    left: Box::new(NodeV2::Leaf(left)),
                    right: Box::new(NodeV2::Leaf(right)),
                    hash: hasher.finalize().into(),
                };
                let node = NodeV2::Branch(branch);
                branches.push(Box::new(node));
                pair.clear();
            }
        }

        // We iteratively pair branches/nodes until we reach depth
        for _ in 1..depth {
            let mut next_branches: Vec<Box<NodeV2<T>>> = vec![];
            let mut pair: Vec<Box<NodeV2<T>>> = vec![];
            for b in branches.into_iter() {
                pair.push(b);
                if pair.len() == 2 {
                    let right = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                    let left = pair.pop().ok_or(anyhow!("unexpected empty pair"))?;
                    let mut hasher = Sha3_256::new();
                    hasher.update(left.hash());
                    hasher.update(right.hash());
                    let branch = BranchV2 {
                        left,
                        right,
                        hash: hasher.finalize().into(),
                    };
                    let node = NodeV2::Branch(branch);
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
            NodeV2::Branch(b) => b,
            NodeV2::Leaf(_) => return Err(anyhow!("unexpected leaf at the root retrieval")),
        };

        Ok(MerkleTreeV2 { depth, root })
    }
}

#[allow(dead_code)]
pub fn verify_proof<T>(
    root: &[u8; 32],
    value: T,
    proof: &[(bool, [u8; 32])],
) -> Result<(), anyhow::Error>
where
    T: AsRef<[u8]>,
{
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

#[cfg(test)]
mod merkletree_tests {
    use super::*;
    use bitvec::prelude::*;

    #[derive(Clone, PartialEq, Debug)]
    struct Stuff {
        v: [u8; 64],
    }

    impl AsRef<[u8]> for Stuff {
        fn as_ref(&self) -> &[u8] {
            &self.v
        }
    }

    impl Default for Stuff {
        fn default() -> Self {
            Stuff { v: [0; 64] }
        }
    }

    #[test]
    fn test_v2() {
        let depth = 10u8;
        let number_of_elements = 2usize.pow(depth.into());
        let mut values = vec![];
        for _ in 0..number_of_elements {
            let inner_value: [u8; 64] = rand::random();
            values.push(Stuff { v: inner_value });
        }
        let tree = MerkleTreeV2::new(depth, &values).unwrap();

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
            let owned_proof: Vec<(bool, [u8; 32])> = selected
                .proof
                .iter()
                .map(|v| (v.0, v.1.to_owned()))
                .collect();
            if let Err(e) = verify_proof(tree.root_hash(), selected.value, &owned_proof) {
                panic!("{e:?}");
            }
        }
    }
}
