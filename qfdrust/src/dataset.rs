use std::collections::{BTreeMap};
use rayon::prelude::*;
use crate::utility;

#[derive(Clone, Debug)]
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
    pub q_value: Option<f64>,
    pub confidence: Option<f64>,
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

    pub fn size(&self) -> usize {
        self.psm_map.len()
    }
    pub fn tdc(&self, method: TDCMethod) -> Vec<PeptideSpectrumMatch> {
        match method {
            TDCMethod::PsmLevel => tdc_psm_level(&self),
            TDCMethod::PeptideLevelPsmOnly => tdc_peptide_only(&self),
            TDCMethod::PeptideLevelPeptideOnly => tdc_peptide_only(&self),
            TDCMethod::PeptideLevelPsmAndPeptide => tdc_psm_and_peptide(&self),
        }
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

/// Get the best target and decoy match for each spectrum and select the best one
///
/// # Arguments
///
/// * `ds` - A PsmDataset
///
/// # Returns
///
/// A vector of PeptideSpectrumMatch
///
fn get_candidates_psm_only(ds: &PsmDataset) -> Vec<&PeptideSpectrumMatch> {

    let mut restult: Vec<&PeptideSpectrumMatch> = Vec::new();

    for spec_id in ds.get_spectra_ids().iter() {
        // for each spectrum, get the best target and decoy match
        let matches = ds.psm_map.get(spec_id).unwrap();
        let maybe_best_target = matches.iter().find(|psm| !psm.decoy);
        let maybe_best_decoy = matches.iter().find(|psm| psm.decoy);

        // if both target and decoy are present, select the best one
        match (maybe_best_target, maybe_best_decoy) {
            (Some(best_target), Some(best_decoy)) => {
                if best_target.score > best_decoy.score {
                    restult.push(best_target);
                }
                // randomly select target or decoy if they have the same score
                if best_target.score == best_decoy.score {
                    if rand::random() {
                        restult.push(best_target);
                    } else {
                        restult.push(best_decoy);
                    }
                }
                if best_target.score < best_decoy.score {
                    restult.push(best_decoy);
                }
            },

            // if only target is present, select it
            (Some(best_target), None) => {
                restult.push(best_target);
            },

            // if only decoy is present, select it
            (None, Some(best_decoy)) => {
                restult.push(best_decoy);
            },

            _ => {},
        }
    }

    restult.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    restult
}

/// Perform target-decoy competition at the PSM level
///
/// # Arguments
///
/// * `ds` - A PsmDataset
///
/// # Returns
///
/// A tuple containing the selected PeptideSpectrumMatch and the q-values, caution,
/// random selection of target or decoy at the same score make this function non-deterministic,
fn tdc_psm_level(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    let candidates = get_candidates_psm_only(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        PeptideSpectrumMatch {
            spec_id: psm.spec_id.clone(),
            peptide_id: psm.peptide_id,
            proteins: psm.proteins.clone(),
            decoy: psm.decoy,
            score: psm.score,
            intensity_ms1: psm.intensity_ms1,
            intensity_ms2: psm.intensity_ms2,
            features: psm.features.clone(),
            q_value: Some(*q),
            confidence: None,
        }
    }).collect()
}

fn tdc_peptide_only(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    todo!("Implement this function")
}

fn tdc_psm_and_peptide(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    todo!("Implement this function")
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