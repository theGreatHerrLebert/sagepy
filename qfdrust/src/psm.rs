use std::collections::BTreeMap;
use rustms::chemistry::formula::calculate_mz;
use rustms::proteomics::peptide::{FragmentType, PeptideProductIonSeriesCollection, PeptideSequence};
use sage_core::ion_series::Kind;
use sage_core::scoring::{Feature};
use serde::{Deserialize, Serialize};
use crate::intensity::FragmentIntensityPrediction;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Psm {
    pub spec_idx: String,
    pub peptide_idx: u32,
    pub proteins: Vec<String>,
    pub hyperscore: f64,
    pub decoy: bool,
    pub sage_feature: Feature,
    pub sequence: Option<PeptideSequence>,
    pub charge: Option<u8>,
    pub mono_mz_calculated: Option<f32>,
    pub mono_mass_observed: Option<f32>,
    pub mono_mass_calculated: Option<f32>,
    pub intensity_ms1: Option<f32>,
    pub intensity_ms2: Option<f32>,
    pub collision_energy: Option<f32>,
    pub collision_energy_calibrated: Option<f32>,
    pub retention_time: Option<f32>,
    pub retention_time_calibrated: Option<f32>,
    pub inverse_ion_mobility: Option<f32>,
    pub inverse_ion_mobility_calibrated: Option<f32>,
    pub prosit_predicted_intensities: Option<Vec<f32>>,
    pub re_score: Option<f64>,
    pub q_value: Option<f64>,
    pub posterior_error_probability: Option<f64>,
    pub external_features: BTreeMap<String, f32>,
}

impl Psm {
    pub fn new(
        spec_idx: String,
        peptide_idx: u32,
        proteins: Vec<String>,
        hyperscore: f64,
        decoy: bool,
        sage_feature: Feature,
        sequence: Option<String>,
        charge: Option<u8>,
        mono_mass_observed: Option<f32>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        collision_energy: Option<f32>,
        collision_energy_calibrated: Option<f32>,
        retention_time: Option<f32>,
        retention_time_calibrated: Option<f32>,
        inverse_ion_mobility: Option<f32>,
        inverse_ion_mobility_calibrated: Option<f32>,
        prosit_predicted_intensities: Option<Vec<f32>>,
        re_score: Option<f64>,
        q_value: Option<f64>,
        posterior_error_probability: Option<f64>,
    ) -> Self {

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
        
        Psm {
            spec_idx,
            peptide_idx,
            proteins,
            hyperscore,
            decoy,
            sage_feature,
            sequence: peptide_sequence,
            charge,
            mono_mz_calculated,
            mono_mass_observed,
            mono_mass_calculated,
            intensity_ms1,
            intensity_ms2,
            collision_energy,
            collision_energy_calibrated,
            retention_time,
            retention_time_calibrated,
            inverse_ion_mobility,
            inverse_ion_mobility_calibrated,
            prosit_predicted_intensities,
            re_score,
            q_value,
            posterior_error_probability,
            external_features: BTreeMap::new(),
        }
    }
    pub fn associate_with_prosit_predicted_intensities(&self, flat_intensities: Vec<f64>) -> Option<PeptideProductIonSeriesCollection> {
        match &self.sequence {
            Some(seq) => Some(seq.associate_with_predicted_intensities(self.charge.unwrap() as i32, FragmentType::B, flat_intensities, false, false)),
            None => None,
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
}