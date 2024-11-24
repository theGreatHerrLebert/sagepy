use std::collections::BTreeMap;
use rustms::chemistry::formula::calculate_mz;
use rustms::proteomics::peptide::{PeptideSequence};
use sage_core::ion_series::Kind;
use sage_core::scoring::{Feature, Fragments};
use serde::{Deserialize, Serialize};
use crate::intensity::{prosit_intensities_to_fragments, FragmentIntensityPrediction};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Psm {
    pub spec_idx: String,
    pub peptide_idx: u32,
    pub proteins: Vec<String>,
    pub sage_feature: Feature,
    pub sequence: Option<PeptideSequence>,
    pub mono_mz_calculated: Option<f32>,
    pub intensity_ms1: Option<f32>,
    pub intensity_ms2: Option<f32>,
    pub collision_energy: Option<f32>,
    pub collision_energy_calibrated: Option<f32>,
    pub retention_time_projected: Option<f32>,
    pub prosit_predicted_intensities: Option<Vec<f32>>,
    pub re_score: Option<f64>,
    pub fragment_intensity_prediction: Option<FragmentIntensityPrediction>,
}

impl Psm {
    pub fn new(
        spec_idx: String,
        peptide_idx: u32,
        proteins: Vec<String>,
        sage_feature: Feature,
        sequence: Option<String>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        collision_energy: Option<f32>,
        collision_energy_calibrated: Option<f32>,
        retention_time_projected: Option<f32>,
        prosit_predicted_intensities: Option<Vec<f32>>,
        re_score: Option<f64>,
    ) -> Self {

        let peptide_sequence = match &sequence {
            Some(seq) => Some(PeptideSequence::new(seq.clone(), Some(peptide_idx as i32))),
            None => None,
        };

        let mono_mz_calculated = match (peptide_sequence.clone(), sage_feature.charge as i32) {
            (Some(seq), charge) => Some(calculate_mz(seq.mono_isotopic_mass(), charge) as f32),
            (_, _) => None,
        };
        
        Psm {
            spec_idx,
            peptide_idx,
            proteins,
            sage_feature,
            sequence: peptide_sequence,
            mono_mz_calculated,
            intensity_ms1,
            intensity_ms2,
            collision_energy,
            collision_energy_calibrated,
            retention_time_projected,
            prosit_predicted_intensities,
            re_score,
            fragment_intensity_prediction: None,
        }
    }

    pub fn set_prosit_intensity_prediction(& mut self, flat_prosit_intensities: Option<Vec<f32>>) {
        self.prosit_predicted_intensities = flat_prosit_intensities.clone();
        match flat_prosit_intensities {
            Some(_intensities) => self.fragment_intensity_prediction = Some(self.get_fragment_intensity_prediction()),
            None => self.fragment_intensity_prediction = None
        }
    }

    pub fn get_fragment_intensity_prediction(&self) -> FragmentIntensityPrediction {
        FragmentIntensityPrediction::new(
            self.sage_feature.fragments.clone().unwrap().intensities,
            self.sage_feature.fragments.clone().unwrap().mz_experimental,
            self.sage_feature.fragments.clone().unwrap().mz_calculated,
            self.sage_feature.fragments.clone().unwrap().charges,
            self.sage_feature.fragments.clone().unwrap().fragment_ordinals,
            self.sage_feature.fragments.clone().unwrap().kinds.iter().map(|kind| *kind == Kind::Y).collect(),
            self.prosit_predicted_intensities.clone().unwrap(),
        )
    }

   pub fn calculate_fragment_intensity_prediction(&mut self) {
        self.fragment_intensity_prediction = Some(self.get_fragment_intensity_prediction());
    }

    pub fn prosit_intensity_to_fragments(&self) -> Option<Fragments> {
        match &self.prosit_predicted_intensities {
            Some(intensities) => Some(prosit_intensities_to_fragments(intensities.clone())),
            None => None,
        }
    }

    pub fn to_dict(&self) -> BTreeMap<String, f64> {
        let mut dict = BTreeMap::new();

        let sage_feature = &self.sage_feature;
        dict.insert("sage_peptide_len".to_string(), sage_feature.peptide_len as f64);
        dict.insert("sage_rank".to_string(), sage_feature.rank as f64);
        dict.insert("sage_expmass".to_string(), sage_feature.expmass as f64);
        dict.insert("sage_calcmass".to_string(), sage_feature.calcmass as f64);
        dict.insert("sage_rt".to_string(), sage_feature.rt as f64);
        dict.insert("sage_aligned_rt".to_string(), sage_feature.aligned_rt as f64);
        dict.insert("sage_predicted_rt".to_string(), sage_feature.predicted_rt as f64);
        dict.insert("sage_delta_rt_model".to_string(), sage_feature.delta_rt_model as f64);
        dict.insert("sage_ims".to_string(), sage_feature.ims as f64);
        dict.insert("sage_predicted_ims".to_string(), sage_feature.predicted_ims as f64);
        dict.insert("sage_delta_ims_model".to_string(), sage_feature.delta_ims_model as f64);
        dict.insert("sage_delta_mass".to_string(), sage_feature.delta_mass as f64);
        dict.insert("sage_isotope_error".to_string(), sage_feature.isotope_error as f64);
        dict.insert("sage_average_ppm".to_string(), sage_feature.average_ppm as f64);
        dict.insert("sage_hyperscore".to_string(), sage_feature.hyperscore);
        dict.insert("sage_delta_next".to_string(), sage_feature.delta_next);
        dict.insert("sage_delta_best".to_string(), sage_feature.delta_best);
        dict.insert("sage_matched_peaks".to_string(), sage_feature.matched_peaks as f64);
        dict.insert("sage_longest_b".to_string(), sage_feature.longest_b as f64);
        dict.insert("sage_longest_y".to_string(), sage_feature.longest_y as f64);
        dict.insert("sage_longest_y_pct".to_string(), sage_feature.longest_y_pct as f64);
        dict.insert("sage_missed_cleavages".to_string(), sage_feature.missed_cleavages as f64);
        dict.insert("sage_matched_intensity_pct".to_string(), sage_feature.matched_intensity_pct as f64);
        dict.insert("sage_scored_candidates".to_string(), sage_feature.scored_candidates as f64);
        dict.insert("sage_poisson".to_string(), sage_feature.poisson);
        dict.insert("sage_discriminant_score".to_string(), sage_feature.discriminant_score as f64);
        dict.insert("sage_posterior_error".to_string(), sage_feature.posterior_error as f64);
        dict.insert("sage_spectrum_q".to_string(), sage_feature.spectrum_q as f64);
        dict.insert("sage_peptide_q".to_string(), sage_feature.peptide_q as f64);
        dict.insert("sage_protein_q".to_string(), sage_feature.protein_q as f64);
        dict.insert("sage_ms2_intensity".to_string(), sage_feature.ms2_intensity as f64);
        dict.insert("decoy".to_string(), self.sage_feature.label as f64);

        dict.insert("intensity_ms1".to_string(), self.intensity_ms1.unwrap_or(0.0) as f64);
        dict.insert("intensity_ms2".to_string(), self.intensity_ms2.unwrap_or(0.0) as f64);
        dict.insert("collision_energy".to_string(), self.collision_energy.unwrap_or(0.0) as f64);
        dict.insert("collision_energy_calibrated".to_string(), self.collision_energy_calibrated.unwrap_or(0.0) as f64);
        dict.insert("retention_time_projected".to_string(), self.retention_time_projected.unwrap_or(0.0) as f64);
        dict.insert("re_score".to_string(), self.re_score.unwrap_or(0.0));

        let intensity_features = self.fragment_intensity_prediction.clone().unwrap();
        dict.insert("intensity_cosine_similarity".to_string(), intensity_features.cosine_similarity(0.0, false).unwrap_or(0.0) as f64);
        dict.insert("intensity_spectral_angle_similarity".to_string(), intensity_features.spectral_angle_similarity(0.0, false) as f64);
        dict.insert("intensity_pearson_correlation".to_string(), intensity_features.pearson_correlation(0.0, false) as f64);
        dict.insert("intensity_spearman_correlation".to_string(), intensity_features.spearman_correlation(0.0, false) as f64);
        dict.insert("intensity_spectral_entropy_similarity".to_string(), intensity_features.spectral_entropy_similarity(0.0, false) as f64);

        dict
    }
}