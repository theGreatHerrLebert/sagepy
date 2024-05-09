use std::collections::{BTreeMap};
use itertools::multizip;
use rustms::chemistry::formula::calculate_mz;
use rustms::proteomics::peptide::{FragmentType, PeptideProductIonSeriesCollection, PeptideSequence};
use serde::{Deserialize, Serialize};
use crate::utility;

#[derive(Clone, Debug)]
pub struct Match {
    pub match_idx: u32,
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

    pub fn from_vectors(spectrum_idx: Vec<String>, match_idx: Vec<u32>, decoy: Vec<bool>, score: Vec<f32>) -> MatchDataset {
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

    pub fn to_vectors(&self) -> (Vec<String>, Vec<u32>, Vec<bool>, Vec<f32>, Vec<f64>) {
        let mut spectrum_idx: Vec<String> = Vec::new();
        let mut match_idx: Vec<u32> = Vec::new();
        let mut decoy: Vec<bool> = Vec::new();
        let mut score: Vec<f32> = Vec::new();
        let mut q_value: Vec<f64> = Vec::new();

        for (spec_idx, matches) in &self.matches {
            for m in matches {
                spectrum_idx.push(spec_idx.clone());
                match_idx.push(m.match_idx);
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PeptideSpectrumMatch {
    pub spec_idx: String,
    pub peptide_idx: u32,
    pub proteins: Vec<String>,
    pub decoy: bool,
    pub hyper_score: f64,
    pub rank: u32,
    pub mono_mz_calculated: Option<f32>,
    pub mono_mass_observed: Option<f32>,
    pub mono_mass_calculated: Option<f32>,
    pub isotope_error: Option<f32>,
    pub average_ppm: Option<f32>,
    pub delta_next: Option<f64>,
    pub delta_best: Option<f64>,
    pub matched_peaks: Option<u32>,
    pub longest_b: Option<u32>,
    pub longest_y: Option<u32>,
    pub longest_y_pct: Option<f32>,
    pub missed_cleavages: Option<u8>,
    pub matched_intensity_pct: Option<f32>,
    pub scored_candidates: Option<u32>,
    pub poisson: Option<f64>,
    pub peptide_sequence: Option<PeptideSequence>,
    pub charge: Option<u8>,
    pub retention_time_observed: Option<f32>,
    pub retention_time_predicted: Option<f32>,
    pub inverse_mobility_observed: Option<f32>,
    pub inverse_mobility_predicted: Option<f32>,
    pub intensity_ms1: Option<f32>,
    pub intensity_ms2: Option<f32>,
    pub q_value: Option<f64>,
    pub collision_energy: Option<f64>,
    pub collision_energy_calibrated: Option<f64>,
    pub re_score: Option<f64>,
    pub cosine_similarity: Option<f32>,
    pub file_name: Option<String>,
}

impl PeptideSpectrumMatch {
    pub fn new(
        spec_idx: String,
        peptide_idx: u32,
        proteins: Vec<String>,
        decoy: bool,
        hyper_score: f64,
        rank: u32,
        mono_mass_observed: Option<f32>,
        isotope_error: Option<f32>,
        average_ppm: Option<f32>,
        delta_next: Option<f64>,
        delta_best: Option<f64>,
        matched_peaks: Option<u32>,
        longest_b: Option<u32>,
        longest_y: Option<u32>,
        longest_y_pct: Option<f32>,
        missed_cleavages: Option<u8>,
        matched_intensity_pct: Option<f32>,
        scored_candidates: Option<u32>,
        poisson: Option<f64>,
        sequence: Option<String>,
        charge: Option<u8>,
        retention_time_observed: Option<f32>,
        retention_time_predicted: Option<f32>,
        inverse_mobility_observed: Option<f32>,
        inverse_mobility_predicted: Option<f32>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        q_value: Option<f64>,
        collision_energy: Option<f64>,
        collision_energy_calibrated: Option<f64>,
        re_score: Option<f64>,
        cosine_similarity: Option<f32>,
        file_name: Option<String>,
    ) -> PeptideSpectrumMatch {

        let peptide_sequence = match &sequence {
            Some(seq) => Some(PeptideSequence::new(seq.clone(), Some(peptide_idx as i32))),
            None => None,
        };

        let mono_mass_calculated = match peptide_sequence.clone() {
            Some(seq) => Some(seq.mono_isotopic_mass() as f32),
            _ => None,
        };

        let mono_mz_calculated = match (peptide_sequence.clone(), charge) {
            (Some(seq), Some(ch)) => Some(calculate_mz(seq.mono_isotopic_mass(), ch as i32) as f32),
            (_, _) => None,
        };

        PeptideSpectrumMatch {
            spec_idx,
            peptide_idx,
            proteins,
            decoy,
            hyper_score,
            rank,
            mono_mz_calculated,
            mono_mass_observed,
            mono_mass_calculated,
            isotope_error,
            average_ppm,
            delta_next,
            delta_best,
            matched_peaks,
            longest_b,
            longest_y,
            longest_y_pct,
            missed_cleavages,
            matched_intensity_pct,
            scored_candidates,
            poisson,
            peptide_sequence,
            charge,
            retention_time_observed,
            retention_time_predicted,
            inverse_mobility_observed,
            inverse_mobility_predicted,
            intensity_ms1,
            intensity_ms2,
            q_value,
            collision_energy,
            collision_energy_calibrated,
            re_score,
            cosine_similarity,
            file_name,
        }
    }
    pub fn associate_with_prosit_predicted_intensities(&self, flat_intensities: Vec<f64>) -> Option<PeptideProductIonSeriesCollection> {
        match &self.peptide_sequence {
            Some(seq) => Some(seq.associate_with_predicted_intensities(self.charge.unwrap() as i32, FragmentType::B, flat_intensities, false, false)),
            None => None,
        }
    }

    pub fn from_json(json: &str) -> PeptideSpectrumMatch {
        serde_json::from_str(json).unwrap()
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
    let mut peptide_map: BTreeMap<(u32, bool), &Match> = BTreeMap::new();

    for psm in candidates {
        let key = (psm.match_idx, psm.decoy);
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

    let mut peptide_map: BTreeMap<u32, (Option<&Match>, Option<&Match>)> = BTreeMap::new();

    for psm in best_targets.iter().chain(best_decoys.iter()) {
        let key = psm.match_idx;
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
    let mut peptide_map: BTreeMap<u32, (Option<&Match>, Option<&Match>)> = BTreeMap::new();

    for psm in best_psms {
        let key = psm.match_idx;
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
            match_idx: psm.match_idx,
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
            match_idx: psm.match_idx,
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
            match_idx: psm.match_idx,
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
            match_idx: psm.match_idx,
            match_identity_candidates: psm.match_identity_candidates.clone(),
            decoy: psm.decoy,
            score: psm.score,
            q_value: Some(*q),
        }
    }).collect()
}

pub fn target_decoy_competition(method: TDCMethod, spectra_idx: Vec<String>, match_idx: Vec<u32>, is_decoy: Vec<bool>, scores: Vec<f32>) -> (Vec<String>, Vec<u32>, Vec<bool>, Vec<f32>, Vec<f64>) {
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