use std::collections::{BTreeMap};
use rayon::prelude::*;

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
    pub fn get_spectra_ids(&self) -> Vec<String> {
        self.psm_map.keys().cloned().collect()
    }

    pub fn get_best_target_psm(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        self.psm_map.get(spec_id).and_then(|psms| psms.iter().find(|psm| !psm.decoy))
    }

    pub fn get_best_decoy_psm(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        self.psm_map.get(spec_id).and_then(|psms| psms.iter().find(|psm| psm.decoy))
    }

    pub fn get_best_match(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        let best_target = self.get_best_target_psm(spec_id);
        let best_decoy = self.get_best_decoy_psm(spec_id);
        // randomly return target or decoy if both are present and have the same score, otherwise return the best one
        match (best_target, best_decoy) {
            (Some(target), Some(decoy)) => {
                if target.score == decoy.score {
                    if rand::random() {
                        Some(target)
                    } else {
                        Some(decoy)
                    }
                } else {
                    Some(target)
                }
            },
            (Some(target), None) => Some(target),
            (None, Some(decoy)) => Some(decoy),
            (None, None) => None,
        }
    }

    pub fn get_best_matches(&self, spec_ids: Vec<String>, num_threads: usize) -> Vec<Option<&PeptideSpectrumMatch>> {
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();

        let result = thread_pool.install(|| {
            spec_ids.into_par_iter().map(|spec_id| self.get_best_match(spec_id.as_str())).collect()
        });

        result
    }

    pub fn size(&self) -> usize {
        self.psm_map.len()
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