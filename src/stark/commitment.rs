use crate::{
    merkletree_v2::{MerkleProof, MerkleTreeV2},
    primefield::PrimeFieldElement,
};

pub fn commit<const N: u64>(
    values: &[Evaluation<N>],
) -> Result<MerkleTreeV2<Evaluation<N>>, anyhow::Error> {
    let mut depth: u8 = 2;
    while 2usize.pow(depth as u32) < values.len() {
        depth += 1;
    }
    MerkleTreeV2::new(depth, values)
}

#[derive(Clone, Debug)]
pub struct Evaluation<const N: u64> {
    pub v: PrimeFieldElement<N>,
    pub rep: [u8; 8],
}

impl<const N: u64> Evaluation<N> {
    pub fn new(v: PrimeFieldElement<N>) -> Self {
        let mut rep = [0u8; 8];
        rep[0..8].clone_from_slice(&v.inner().to_le_bytes());
        Self { v, rep }
    }
}

impl<const N: u64> AsRef<[u8]> for Evaluation<N> {
    fn as_ref(&self) -> &[u8] {
        &self.rep
    }
}

impl<const N: u64> Default for Evaluation<N> {
    fn default() -> Self {
        Self {
            v: 0.into(),
            rep: [0u8; 8],
        }
    }
}

impl<const N: u64> MerkleTreeV2<Evaluation<N>> {
    pub fn select_commitment(&self, unit_index: usize) -> Result<Commitment<N>, anyhow::Error> {
        let selected_value = self.get(&self.index_to_selector(unit_index))?;
        Ok(Commitment {
            unit_index,
            evaluation: selected_value.value.clone(),
            proof: selected_value.proof,
        })
    }
}

pub struct Commitment<const N: u64> {
    pub unit_index: usize,
    pub evaluation: Evaluation<N>,
    pub proof: MerkleProof,
}
