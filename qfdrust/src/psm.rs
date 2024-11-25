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

    pub fn get_feature_vector(&self) -> Vec<f64> {

        let sage_feature = &self.sage_feature;
        let mut feature_vector = Vec::new();
        feature_vector.push(sage_feature.expmass as f64);
        feature_vector.push(sage_feature.calcmass as f64);
        feature_vector.push(sage_feature.charge as f64);
        feature_vector.push(sage_feature.rt as f64);
        feature_vector.push(sage_feature.aligned_rt as f64);
        feature_vector.push(sage_feature.predicted_rt as f64);
        feature_vector.push(sage_feature.delta_rt_model as f64);
        feature_vector.push(sage_feature.ims as f64);
        feature_vector.push(sage_feature.predicted_ims as f64);
        feature_vector.push(sage_feature.delta_ims_model as f64);
        feature_vector.push(sage_feature.delta_mass as f64);
        feature_vector.push(sage_feature.isotope_error as f64);
        feature_vector.push(sage_feature.average_ppm as f64);
        feature_vector.push(sage_feature.hyperscore);
        feature_vector.push(sage_feature.delta_next);
        feature_vector.push(sage_feature.delta_best);
        feature_vector.push(sage_feature.matched_peaks as f64);
        feature_vector.push(sage_feature.longest_b as f64);
        feature_vector.push(sage_feature.longest_y as f64);
        feature_vector.push(sage_feature.longest_y_pct as f64);
        feature_vector.push(sage_feature.missed_cleavages as f64);
        feature_vector.push(sage_feature.matched_intensity_pct as f64);
        feature_vector.push(sage_feature.scored_candidates as f64);
        feature_vector.push(sage_feature.poisson);
        feature_vector.push(sage_feature.discriminant_score as f64);
        feature_vector.push(sage_feature.posterior_error as f64);
        feature_vector.push(sage_feature.ms2_intensity as f64);
        feature_vector.push(sage_feature.rank as f64);

        feature_vector.push(self.intensity_ms1.unwrap_or(0.0) as f64);
        feature_vector.push(self.intensity_ms2.unwrap_or(0.0) as f64);
        feature_vector.push(self.collision_energy.unwrap_or(0.0) as f64);
        feature_vector.push(self.collision_energy_calibrated.unwrap_or(0.0) as f64);
        feature_vector.push(self.retention_time_projected.unwrap_or(0.0) as f64);

        let intensity_features = self.fragment_intensity_prediction.clone();

        match intensity_features {
            Some(intensity_features) => {
                let features = intensity_features.get_feature_vector(0.00001, false);
                for feature in features {
                    feature_vector.push(feature as f64);
                }
            },

            None => {
                for _ in 0..5 {
                    feature_vector.push(0.0);
                }
            }
        }

        feature_vector.push(sage_feature.label as f64);
        feature_vector.push(sage_feature.spectrum_q as f64);
        feature_vector.push(sage_feature.peptide_q as f64);
        feature_vector.push(sage_feature.protein_q as f64);

        feature_vector
    }

    pub fn get_feature_names(&self) -> Vec<&str> {
        vec![
            "expmass",
            "calcmass",
            "charge",
            "rt",
            "aligned_rt",
            "predicted_rt",
            "delta_rt_model",
            "ims",
            "predicted_ims",
            "delta_ims_model",
            "delta_mass",
            "isotope_error",
            "average_ppm",
            "hyperscore",
            "delta_next",
            "delta_best",
            "matched_peaks",
            "longest_b",
            "longest_y",
            "longest_y_pct",
            "missed_cleavages",
            "matched_intensity_pct",
            "scored_candidates",
            "poisson",
            "discriminant_score",
            "posterior_error",
            "ms2_intensity",
            "rank",
            "intensity_ms1",
            "intensity_ms2",
            "collision_energy",
            "collision_energy_calibrated",
            "retention_time_projected",
            "cosine_similarity",
            "spectral_angle_similarity",
            "pearson_correlation",
            "spearman_correlation",
            "spectral_entropy_similarity",
            "decoy",
            "spectrum_q",
            "peptide_q",
            "protein_q",
        ]
    }
}