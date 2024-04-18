use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct PeptideSpectrumMatch {
    pub spec_id: String,
    pub peptide_id: u32,
    pub proteins: Vec<String>,
    pub decoy: bool,
    pub score: f64,
}

pub struct PsmDataset {
    pub psm_map: HashMap<String, Vec<PeptideSpectrumMatch>>,
}