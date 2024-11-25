use std::collections::{BTreeMap};
use itertools::multizip;
use crate::utility;

#[derive(Clone, Debug)]
pub struct Match {
    pub match_idx: String,
    pub spectrum_idx: String,
    pub match_identity_candidates: Option<Vec<String>>,
    pub decoy: bool,
    pub score: f32,
    pub q_value: Option<f64>,
}

pub struct MatchDataset {
    pub matches: BTreeMap<String, Vec<Match>>,
}

impl MatchDataset {
    pub fn new(mut map: BTreeMap<String, Vec<Match>>) -> MatchDataset {
        // score should be descending
        for (_, matches) in &mut map {
            matches.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        }

        map.retain(|_, v| !v.is_empty());

        MatchDataset {
            matches: map,
        }
    }

    pub fn get_spectra_ids(&self) -> Vec<String> {
        self.matches.keys().cloned().collect()
    }

    pub fn from_collection(collection: Vec<Match>) -> MatchDataset {
        let mut map: BTreeMap<String, Vec<Match>> = BTreeMap::new();

        for m in collection {
            let entry = map.entry(m.spectrum_idx.clone()).or_insert(Vec::new());
            entry.push(m);
        }

        MatchDataset::new(map)
    }

    pub fn from_vectors(spectrum_idx: Vec<String>, match_idx: Vec<String>, decoy: Vec<bool>, score: Vec<f32>) -> MatchDataset {
        let mut map: BTreeMap<String, Vec<Match>> = BTreeMap::new();
        for (spec_idx, match_idx, score, d) in multizip((spectrum_idx, match_idx, score, decoy)) {
            let entry = map.entry(spec_idx.clone()).or_insert(Vec::new());
            entry.push(Match {
                score,
                match_idx,
                spectrum_idx: spec_idx,
                decoy: d,
                match_identity_candidates: None,
                q_value: None,
            });
        }
        MatchDataset::new(map)
    }

    pub fn to_vectors(&self) -> (Vec<String>, Vec<String>, Vec<bool>, Vec<f32>, Vec<f64>) {
        let mut spectrum_idx: Vec<String> = Vec::new();
        let mut match_idx: Vec<String> = Vec::new();
        let mut decoy: Vec<bool> = Vec::new();
        let mut score: Vec<f32> = Vec::new();
        let mut q_value: Vec<f64> = Vec::new();

        for (spec_idx, matches) in &self.matches {
            for m in matches {
                spectrum_idx.push(spec_idx.clone());
                match_idx.push(m.match_idx.clone());
                decoy.push(m.decoy);
                score.push(m.score);
                q_value.push(m.q_value.unwrap_or(1.0));
            }
        }

        (spectrum_idx, match_idx, decoy, score, q_value)
    }

    pub fn get_best_target_match(&self, spec_id: &str) -> Option<&Match> {
        self.matches.get(spec_id).unwrap().iter().find(|m| !m.decoy)
    }

    pub fn get_best_decoy_match(&self, spec_id: &str) -> Option<&Match> {
        self.matches.get(spec_id).unwrap().iter().find(|m| m.decoy)
    }

    pub fn get_best_target_matches(&self) -> Vec<&Match> {
        let maybe_matches: Vec<_> = self.get_spectra_ids().iter().map(|spec_id| self.get_best_target_match(spec_id)).collect();
        maybe_matches.iter().filter_map(|m| *m).collect()
    }

    pub fn get_best_decoy_matches(&self) -> Vec<&Match> {
        let maybe_matches: Vec<_> = self.get_spectra_ids().iter().map(|spec_id| self.get_best_decoy_match(spec_id)).collect();
        maybe_matches.iter().filter_map(|m| *m).collect()
    }

    pub fn get_best_match(&self, spec_id: &str) -> Option<&Match> {
        let matches = self.matches.get(spec_id).unwrap();
        let maybe_best_target = matches.iter().find(|m| !m.decoy);
        let maybe_best_decoy = matches.iter().find(|m| m.decoy);

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
    pub fn get_best_matches(&self) -> Vec<&Match> {
        self.get_spectra_ids().iter().filter_map(|spec_id| self.get_best_match(spec_id)).collect()
    }
}

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

fn get_candidate_psm_match(ds: &MatchDataset) -> Vec<&Match> {
    let mut result: Vec<&Match> = Vec::new();

    for spec_id in ds.get_spectra_ids().iter() {
        let matches = ds.matches.get(spec_id).unwrap();
        let maybe_best_target = matches.iter().find(|m| !m.decoy);
        let maybe_best_decoy = matches.iter().find(|m| m.decoy);

        match (maybe_best_target, maybe_best_decoy) {
            (Some(best_target), Some(best_decoy)) => {
                if best_target.score > best_decoy.score {
                    result.push(best_target);
                }
                if best_target.score == best_decoy.score {
                    if rand::random() {
                        result.push(best_target);
                    } else {
                        result.push(best_decoy);
                    }
                }
                if best_target.score < best_decoy.score {
                    result.push(best_decoy);
                }
            },
            (Some(best_target), None) => {
                result.push(best_target);
            },
            (None, Some(best_decoy)) => {
                result.push(best_decoy);
            },
            _ => {},
        }
    }

    result.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    result
}

fn get_candidates_peptide_psm_only_match(ds: &MatchDataset) -> Vec<&Match> {
    let candidates = get_candidate_psm_match(ds);
    let mut peptide_map: BTreeMap<(String, bool), &Match> = BTreeMap::new();

    for psm in candidates {
        let key = (psm.match_idx.clone(), psm.decoy);
        let entry = peptide_map.entry(key);
        let best_psm = entry.or_insert(psm);
        if psm.score > best_psm.score || (psm.score == best_psm.score && rand::random()) {
            *best_psm = psm;
        }
    }

    let mut result: Vec<&Match> = Vec::new();

    for (_, psm) in peptide_map {
        result.push(psm);
    }

    result
}

fn get_candidates_peptide_peptide_only_match(ds: &MatchDataset) -> Vec<&Match> {
    let best_targets = ds.get_best_target_matches();
    let best_decoys = ds.get_best_decoy_matches();

    let mut peptide_map: BTreeMap<String, (Option<&Match>, Option<&Match>)> = BTreeMap::new();

    for psm in best_targets.iter().chain(best_decoys.iter()) {
        let key = psm.match_idx.clone();
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

    let mut result: Vec<&Match> = Vec::new();

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

fn get_candidates_peptide_psm_peptide_match(ds: &MatchDataset) -> Vec<&Match> {
    let best_psms = ds.get_best_matches();
    let mut peptide_map: BTreeMap<String, (Option<&Match>, Option<&Match>)> = BTreeMap::new();

    for psm in best_psms {
        let key = psm.match_idx.clone();
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

    let mut result: Vec<&Match> = Vec::new();

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

fn tdc_psm_match(ds: &MatchDataset) -> Vec<Match> {
    let candidates = get_candidate_psm_match(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score as f64).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);
    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        Match {
            spectrum_idx: psm.spectrum_idx.clone(),
            match_idx: psm.match_idx.clone(),
            match_identity_candidates: psm.match_identity_candidates.clone(),
            decoy: psm.decoy,
            score: psm.score,
            q_value: Some(*q),
        }
    }).collect()
}

fn tdc_peptide_psm_only_match(ds: &MatchDataset) -> Vec<Match> {
    let candidates = get_candidates_peptide_psm_only_match(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score as f64).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);
    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        Match {
            spectrum_idx: psm.spectrum_idx.clone(),
            match_idx: psm.match_idx.clone(),
            match_identity_candidates: psm.match_identity_candidates.clone(),
            decoy: psm.decoy,
            score: psm.score,
            q_value: Some(*q),
        }
    }).collect()
}

fn tdc_peptide_peptide_only_match(ds: &MatchDataset) -> Vec<Match> {
    let candidates = get_candidates_peptide_peptide_only_match(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score as f64).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);
    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        Match {
            spectrum_idx: psm.spectrum_idx.clone(),
            match_idx: psm.match_idx.clone(),
            match_identity_candidates: psm.match_identity_candidates.clone(),
            decoy: psm.decoy,
            score: psm.score,
            q_value: Some(*q),
        }
    }).collect()
}
fn tdc_peptide_psm_peptide_match(ds: &MatchDataset) -> Vec<Match> {
    let candidates = get_candidates_peptide_psm_peptide_match(ds);
    let scores: Vec<f64> = candidates.iter().map(|psm| psm.score as f64).collect();
    let targets: Vec<bool> = candidates.iter().map(|psm| !psm.decoy).collect();
    let q_values: Vec<f64> = utility::target_decoy_competition(&scores, &targets, true);

    candidates.iter().zip(q_values.iter()).map(|(psm, q)| {
        Match {
            spectrum_idx: psm.spectrum_idx.clone(),
            match_idx: psm.match_idx.clone(),
            match_identity_candidates: psm.match_identity_candidates.clone(),
            decoy: psm.decoy,
            score: psm.score,
            q_value: Some(*q),
        }
    }).collect()
}

pub fn target_decoy_competition(method: TDCMethod, spectra_idx: Vec<String>, match_idx: Vec<String>, is_decoy: Vec<bool>, scores: Vec<f32>) -> (Vec<String>, Vec<String>, Vec<bool>, Vec<f32>, Vec<f64>) {
    let ds = MatchDataset::from_vectors(spectra_idx, match_idx, is_decoy, scores);

    let result = match method {
        TDCMethod::PsmLevel => tdc_psm_match(&ds),
        TDCMethod::PeptideLevelPsmOnly => tdc_peptide_psm_only_match(&ds),
        TDCMethod::PeptideLevelPeptideOnly => tdc_peptide_peptide_only_match(&ds),
        TDCMethod::PeptideLevelPsmPeptide => tdc_peptide_psm_peptide_match(&ds),
    };

    let (spectrum_idx, match_idx, decoy, score, q_value) = MatchDataset::from_collection(result).to_vectors();
    (spectrum_idx, match_idx, decoy, score, q_value)
}