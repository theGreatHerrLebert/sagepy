use std::collections::{BTreeMap};

#[derive(Clone, Debug)]
pub struct PeptideSpectrumMatch {
    pub spec_id: String,
    pub peptide_id: u32,
    pub proteins: Vec<String>,
    pub decoy: bool,
    pub score: f64,
    pub intensity_ms1: Option<f64>,
    pub intensity_ms2: Option<f64>,
    pub features: Option<Vec<(String, f64)>>,
}

pub struct PsmDataset {
    pub psm_map: BTreeMap<String, Vec<PeptideSpectrumMatch>>,
}

impl PsmDataset {
    pub fn new(mut map: BTreeMap<String, Vec<PeptideSpectrumMatch>>) -> PsmDataset {
        for (_, psms) in &mut map {
            psms.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
        }
        PsmDataset {
            psm_map: map,
        }
    }
    pub fn size(&self) -> usize {
        self.psm_map.len()
    }
    pub fn inverted_index(&self) -> BTreeMap<(u32, bool), Vec<&PeptideSpectrumMatch>> {
        let mut inverted_index: BTreeMap<(u32, bool), Vec<&PeptideSpectrumMatch>> = BTreeMap::new();

        for (_, psms) in &self.psm_map {
            for psm in psms {
                inverted_index.entry((psm.peptide_id, psm.decoy)).or_insert(Vec::new()).push(psm);
            }
        }

        inverted_index
    }
}