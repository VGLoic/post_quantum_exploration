use crate::{
    merkletree_v2::{MerkleProof, MerkleTreeV2},
    primefield::PrimeFieldElement,
};

pub fn build_tree<const N: u32>(
    values: &[Evaluation<N>],
) -> Result<MerkleTreeV2<Evaluation<N>>, anyhow::Error> {
    let mut depth: u8 = 2;
    while 2usize.pow(depth as u32) < values.len() {
        depth += 1;
    }
    MerkleTreeV2::new(depth, values)
}

#[derive(Clone, Debug)]
pub struct Evaluation<const N: u32> {
    pub v: PrimeFieldElement<N>,
    pub rep: [u8; 4],
}

impl<const N: u32> Evaluation<N> {
    pub fn new(v: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 4];
        rep[0..4].clone_from_slice(&v.inner().to_le_bytes());
        Self { v, rep }
    }
}

impl<const N: u32> AsRef<[u8]> for Evaluation<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u32> Default for Evaluation<N> {
    fn default() -> Self {
        Self {
            v: 0.into(),
            rep: [0u8; 4],
        }
    }
}

impl<const N: u32> MerkleTreeV2<Evaluation<N>> {
    pub fn select_commitment(&self, unit_index: usize) -> Result<Commitment<N>, anyhow::Error> {
        let selected_value = self.get(&self.index_to_selector(unit_index))?;
        Ok(Commitment {
            unit_index,
            evaluation: selected_value.value.clone(),
            proof: selected_value.proof,
        })
    }
}

pub struct Commitment<const N: u32> {
    pub unit_index: usize,
    pub evaluation: Evaluation<N>,
    pub proof: MerkleProof,
}
