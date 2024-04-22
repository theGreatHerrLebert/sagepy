use std::collections::{BTreeMap};
use crate::utility;

#[derive(Clone, Debug)]
pub enum TDCMethod {
    PsmLevel,
    PeptideLevelPsmOnly,
    PeptideLevelPeptideOnly,
    PeptideLevelPsmPeptide,
}

impl TDCMethod {
    pub fn from_str(s: &str) -> TDCMethod {
        match s {
            "psm" => TDCMethod::PsmLevel,
            "peptide_psm_only" => TDCMethod::PeptideLevelPsmOnly,
            "peptide_peptide_only" => TDCMethod::PeptideLevelPeptideOnly,
            "peptide_psm_peptide" => TDCMethod::PeptideLevelPsmPeptide,
            _ => panic!("Invalid TDC method"),
        }
    }

    pub fn from_int(i: u32) -> TDCMethod {
        match i {
            0 => TDCMethod::PsmLevel,
            1 => TDCMethod::PeptideLevelPsmOnly,
            2 => TDCMethod::PeptideLevelPeptideOnly,
            3 => TDCMethod::PeptideLevelPsmPeptide,
            _ => panic!("Invalid TDC method"),
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            TDCMethod::PsmLevel => "PsmLevel",
            TDCMethod::PeptideLevelPsmOnly => "PeptideLevelPsmOnly",
            TDCMethod::PeptideLevelPeptideOnly => "PeptideLevelPeptideOnly",
            TDCMethod::PeptideLevelPsmPeptide => "PeptideLevelPsmAndPeptide",
        }
    }
}

#[derive(Clone, Debug)]
pub struct PeptideSpectrumMatch {
    pub spec_id: String,
    pub peptide_id: u32,
    pub sequence: String,
    pub proteins: Vec<String>,
    pub decoy: bool,
    pub score: f64,
    pub intensity_ms1: Option<f64>,
    pub intensity_ms2: Option<f64>,
    pub features: Option<Vec<(String, f64)>>,
    pub q_value: Option<f64>,
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

    pub fn get_candidates(&self, method: TDCMethod) -> Vec<&PeptideSpectrumMatch> {
        match method {
            TDCMethod::PsmLevel => get_candidates_psm(&self),
            TDCMethod::PeptideLevelPsmOnly => get_candidates_peptide_psm_only(&self),
            TDCMethod::PeptideLevelPeptideOnly => get_candidates_peptide_peptide_only(&self),
            TDCMethod::PeptideLevelPsmPeptide => get_candidates_peptide_psm_peptide(&self),
        }
    }

    pub fn tdc(&self, method: TDCMethod) -> Vec<PeptideSpectrumMatch> {

        let mut result = match method {
            TDCMethod::PsmLevel => tdc_psm(&self),
            TDCMethod::PeptideLevelPsmOnly => tdc_peptide_psm_only(&self),
            TDCMethod::PeptideLevelPeptideOnly => tdc_peptide_peptide_only(&self),
            TDCMethod::PeptideLevelPsmPeptide => tdc_peptide_psm_peptide(&self),
        };

       // sort by q value ascending
        result.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap());

        result
    }
    pub fn get_best_target_psm(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        self.psm_map.get(spec_id).unwrap().iter().find(|psm| !psm.decoy)
    }
    pub fn get_best_decoy_psm(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        self.psm_map.get(spec_id).unwrap().iter().find(|psm| psm.decoy)
    }

    pub fn get_best_target_psms(&self) -> Vec<&PeptideSpectrumMatch> {
        let maybe_matches: Vec<_> = self.get_spectra_ids().iter().map(|spec_id| self.get_best_target_psm(spec_id)).collect();
        maybe_matches.iter().filter_map(|m| *m).collect()
    }

    pub fn get_best_decoy_psms(&self) -> Vec<&PeptideSpectrumMatch> {
        let maybe_matches: Vec<_> = self.get_spectra_ids().iter().map(|spec_id| self.get_best_decoy_psm(spec_id)).collect();
        maybe_matches.iter().filter_map(|m| *m).collect()
    }

    pub fn get_best_psm(&self, spec_id: &str) -> Option<&PeptideSpectrumMatch> {
        let matches = self.psm_map.get(spec_id).unwrap();
        let maybe_best_target = matches.iter().find(|psm| !psm.decoy);
        let maybe_best_decoy = matches.iter().find(|psm| psm.decoy);

        match (maybe_best_target, maybe_best_decoy) {
            (Some(best_target), Some(best_decoy)) => {
                if best_target.score > best_decoy.score {
                    Some(best_target)
                }

                else if best_target.score == best_decoy.score {
                    if rand::random() {
                        Some(best_target)
                    } else {
                        Some(best_decoy)
                    }
                }

                else {
                    Some(best_decoy)
                }
            },
            (Some(best_target), None) => Some(best_target),
            (None, Some(best_decoy)) => Some(best_decoy),
            _ => None,
        }
    }
    pub fn get_best_psms(&self) -> Vec<&PeptideSpectrumMatch> {
        self.get_spectra_ids().iter().filter_map(|spec_id| self.get_best_psm(spec_id)).collect()
    }

    pub fn inverted_index(&self) -> BTreeMap<(u32, bool), Vec<&PeptideSpectrumMatch>> {
        let mut peptide_map: BTreeMap<(u32, bool), Vec<&PeptideSpectrumMatch>> = BTreeMap::new();
        for psms in self.psm_map.values() {
            for psm in psms {
                let key = (psm.peptide_id, psm.decoy);
                let entry = peptide_map.entry(key).or_insert(Vec::new());
                entry.push(psm);
            }
        }
        peptide_map
    }
}

// TODO: clean and remove duplicated code

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
fn get_candidates_psm(ds: &PsmDataset) -> Vec<&PeptideSpectrumMatch> {

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

fn get_candidates_peptide_psm_only(ds: &PsmDataset) -> Vec<&PeptideSpectrumMatch> {

    let candidates = get_candidates_psm(ds);
    let mut peptide_map: BTreeMap<(u32, bool), &PeptideSpectrumMatch> = BTreeMap::new();

    // for each peptide, select the best target and decoy match
    for psm in candidates {
        let key = (psm.peptide_id, psm.decoy);
        let entry = peptide_map.entry(key);
        let best_psm = entry.or_insert(psm);
        if psm.score > best_psm.score || (psm.score == best_psm.score && rand::random()) {
            *best_psm = psm;
        }
    }

    let mut result: Vec<&PeptideSpectrumMatch> = Vec::new();

    for (_, psm) in peptide_map {
        result.push(psm);
    }

    result
}

fn get_candidates_peptide_peptide_only(ds: &PsmDataset) -> Vec<&PeptideSpectrumMatch> {

    let best_targets = ds.get_best_target_psms();
    let best_decoys = ds.get_best_decoy_psms();

    let mut peptide_map: BTreeMap<u32, (Option<&PeptideSpectrumMatch>, Option<&PeptideSpectrumMatch>)> = BTreeMap::new();

    for psm in best_targets.iter().chain(best_decoys.iter()) {
        let key = psm.peptide_id;
        let entry = peptide_map.entry(key);
        let (best_target, best_decoy) = entry.or_insert((None, None));
        if psm.decoy {
            if best_decoy.is_none() || psm.score > best_decoy.unwrap().score {
                *best_decoy = Some(psm);
            }
        } else {
            if best_target.is_none() || psm.score > best_target.unwrap().score {
                *best_target = Some(psm);
            }
        }
    }

    let mut result: Vec<&PeptideSpectrumMatch> = Vec::new();

    for (_, (best_target, best_decoy)) in peptide_map {
        match (best_target, best_decoy) {
            (Some(target), Some(decoy)) => {
                if target.score > decoy.score {
                    result.push(target);
                }
                if target.score == decoy.score {
                    if rand::random() {
                        result.push(target);
                    } else {
                        result.push(decoy);
                    }
                }
                else {
                    result.push(decoy);
                }
            },
            (Some(target), None) => {
                result.push(target);
            },
            (None, Some(decoy)) => {
                result.push(decoy);
            },
            _ => {},
        }
    }
    result
}

fn get_candidates_peptide_psm_peptide(ds: &PsmDataset) -> Vec<&PeptideSpectrumMatch> {
    let best_psms = ds.get_best_psms();
    let mut peptide_map: BTreeMap<u32, (Option<&PeptideSpectrumMatch>, Option<&PeptideSpectrumMatch>)> = BTreeMap::new();

    for psm in best_psms {
        let key = psm.peptide_id;
        let entry = peptide_map.entry(key);
        let (best_target, best_decoy) = entry.or_insert((None, None));
        if psm.decoy {
            if best_decoy.is_none() || psm.score > best_decoy.unwrap().score {
                *best_decoy = Some(psm);
            }
        } else {
            if best_target.is_none() || psm.score > best_target.unwrap().score {
                *best_target = Some(psm);
            }
        }
    }

    let mut result: Vec<&PeptideSpectrumMatch> = Vec::new();

    for (_, (best_target, best_decoy)) in peptide_map {
        match (best_target, best_decoy) {
            (Some(target), Some(decoy)) => {
                if target.score > decoy.score {
                    result.push(target);
                }
                if target.score == decoy.score {
                    if rand::random() {
                        result.push(target);
                    } else {
                        result.push(decoy);
                    }
                }
                else {
                    result.push(decoy);
                }
            },
            (Some(target), None) => {
                result.push(target);
            },
            (None, Some(decoy)) => {
                result.push(decoy);
            },
            _ => {},
        }
    }
    result
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
fn tdc_psm(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    let candidates = get_candidates_psm(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        PeptideSpectrumMatch {
            spec_id: psm.spec_id.clone(),
            peptide_id: psm.peptide_id,
            sequence: psm.sequence.clone(),
            proteins: psm.proteins.clone(),
            decoy: psm.decoy,
            score: psm.score,
            intensity_ms1: psm.intensity_ms1,
            intensity_ms2: psm.intensity_ms2,
            features: psm.features.clone(),
            q_value: Some(*q),
        }
    }).collect()
}

fn tdc_peptide_psm_only(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {

    let candidates = get_candidates_peptide_psm_only(ds);

    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        PeptideSpectrumMatch {
            spec_id: psm.spec_id.clone(),
            peptide_id: psm.peptide_id,
            sequence: psm.sequence.clone(),
            proteins: psm.proteins.clone(),
            decoy: psm.decoy,
            score: psm.score,
            intensity_ms1: psm.intensity_ms1,
            intensity_ms2: psm.intensity_ms2,
            features: psm.features.clone(),
            q_value: Some(*q),
        }
    }).collect()
}

fn tdc_peptide_peptide_only(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    let candidates = get_candidates_peptide_peptide_only(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        PeptideSpectrumMatch {
            spec_id: psm.spec_id.clone(),
            peptide_id: psm.peptide_id,
            sequence: psm.sequence.clone(),
            proteins: psm.proteins.clone(),
            decoy: psm.decoy,
            score: psm.score,
            intensity_ms1: psm.intensity_ms1,
            intensity_ms2: psm.intensity_ms2,
            features: psm.features.clone(),
            q_value: Some(*q),
        }
    }).collect()
}

fn tdc_peptide_psm_peptide(ds: &PsmDataset) -> Vec<PeptideSpectrumMatch> {
    let candidates = get_candidates_peptide_psm_peptide(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        PeptideSpectrumMatch {
            spec_id: psm.spec_id.clone(),
            peptide_id: psm.peptide_id,
            sequence: psm.sequence.clone(),
            proteins: psm.proteins.clone(),
            decoy: psm.decoy,
            score: psm.score,
            intensity_ms1: psm.intensity_ms1,
            intensity_ms2: psm.intensity_ms2,
            features: psm.features.clone(),
            q_value: Some(*q),
        }
    }).collect()
}