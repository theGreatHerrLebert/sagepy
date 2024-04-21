use std::collections::{BTreeMap};

pub enum TDCMethod {
    PsmLevel,
    PeptideLevelPsmOnly,
    PeptideLevelPeptideOnly,
    PeptideLevelPsmAndPeptide,
}

impl TDCMethod {
    pub fn from_str(s: &str) -> TDCMethod {
        match s {
            "psm" => TDCMethod::PsmLevel,
            "peptide_psm_only" => TDCMethod::PeptideLevelPsmOnly,
            "peptide_peptide_only" => TDCMethod::PeptideLevelPeptideOnly,
            "peptide_psm_and_peptide" => TDCMethod::PeptideLevelPsmAndPeptide,
            _ => panic!("Invalid TDC method"),
        }
    }

    pub fn from_int(i: u32) -> TDCMethod {
        match i {
            0 => TDCMethod::PsmLevel,
            1 => TDCMethod::PeptideLevelPsmOnly,
            2 => TDCMethod::PeptideLevelPeptideOnly,
            3 => TDCMethod::PeptideLevelPsmAndPeptide,
            _ => panic!("Invalid TDC method"),
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            TDCMethod::PsmLevel => "PsmLevel",
            TDCMethod::PeptideLevelPsmOnly => "PeptideLevelPsmOnly",
            TDCMethod::PeptideLevelPeptideOnly => "PeptideLevelPeptideOnly",
            TDCMethod::PeptideLevelPsmAndPeptide => "PeptideLevelPsmAndPeptide",
        }
    }
}

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
        // score should be descending
        for (_, psms) in &mut map {
            psms.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        }

        map.retain(|_, v| !v.is_empty());

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

    pub fn assign_confidence(&self, method: TDCMethod, threshold: f64, eval_fdr: f64) -> Vec<f64> {
        match method {
            TDCMethod::PsmLevel => assign_confidence_psm_level(&self, threshold, eval_fdr),
            TDCMethod::PeptideLevelPsmOnly => assign_confidence_peptide_level_psm_only(&self, threshold, eval_fdr),
            TDCMethod::PeptideLevelPeptideOnly => assign_confidence_peptide_level_peptide_only(&self, threshold, eval_fdr),
            TDCMethod::PeptideLevelPsmAndPeptide => assign_confidence_peptide_level_psm_and_peptide(&self, threshold, eval_fdr),
        }
    }
}

fn assign_confidence_psm_level(psm_dataset: &PsmDataset, threshold: f64, eval_fdr: f64) -> Vec<f64> {
    assert!(threshold >= 0.0 && threshold <= 1.0, "Threshold should be between 0 and 1");
    todo!("Implement this function")
}

fn assign_confidence_peptide_level_psm_only(psm_dataset: &PsmDataset, threshold: f64, eval_fdr: f64) -> Vec<f64> {
    assert!(threshold >= 0.0 && threshold <= 1.0, "Threshold should be between 0 and 1");
    todo!("Implement this function")
}

fn assign_confidence_peptide_level_peptide_only(psm_dataset: &PsmDataset, threshold: f64, eval_fdr: f64) -> Vec<f64> {
    assert!(threshold >= 0.0 && threshold <= 1.0, "Threshold should be between 0 and 1");
    todo!("Implement this function")
}

fn assign_confidence_peptide_level_psm_and_peptide(psm_dataset: &PsmDataset, threshold: f64, eval_fdr: f64) -> Vec<f64> {
    assert!(threshold >= 0.0 && threshold <= 1.0, "Threshold should be between 0 and 1");
    todo!("Implement this function")
}