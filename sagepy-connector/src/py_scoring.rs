use std::collections::{BTreeMap, HashSet};
use itertools::Itertools;
use pyo3::prelude::*;
use qfdrust::dataset::{PeptideSpectrumMatch};
use qfdrust::intensity::FragmentIntensityPrediction;
use qfdrust::psm::Psm;
use crate::utilities::sage_sequence_to_unimod_sequence;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use sage_core::ion_series::Kind;

use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_mass::PyTolerance;
use crate::py_spectrum::{PyProcessedSpectrum};
use sage_core::scoring::{Feature, Scorer, Fragments, ScoreType};
use sage_core::scoring::ScoreType::{OpenMSHyperScore, SageHyperScore};
use serde::{Deserialize, Serialize};
use crate::py_intensity::PyFragmentIntensityPrediction;
use crate::py_ion_series::PyKind;
use crate::py_utility::{flat_prosit_array_to_fragments_map, py_fragments_to_fragments_map};

#[pyclass]
#[derive(Clone, Serialize)]
pub struct PyPsm {
    pub inner: Psm,
}

#[pymethods]
impl PyPsm {
    #[new]
    pub fn new(
        spec_idx: String,
        peptide_idx: u32,
        proteins: Vec<String>,
        hyperscore: f64,
        decoy: bool,
        sage_feature: PyFeature,
        sequence: Option<String>,
        charge: Option<u8>,
        mono_mass_observed: Option<f32>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        collision_energy: Option<f32>,
        collision_energy_calibrated: Option<f32>,
        retention_time: Option<f32>,
        retention_time_calibrated: Option<f32>,
        retention_time_projected: Option<f32>,
        inverse_ion_mobility: Option<f32>,
        inverse_ion_mobility_calibrated: Option<f32>,
        prosit_predicted_intensities: Option<Vec<f32>>,
        re_score: Option<f64>,
        q_value: Option<f64>,
        posterior_error_probability: Option<f64>,
        py_fragment_intensity_prediction: Option<PyFragmentIntensityPrediction>,
    ) -> Self {
        PyPsm {
            inner: Psm::new(
                spec_idx,
                peptide_idx,
                proteins,
                hyperscore,
                decoy,
                sage_feature.inner.clone(),
                sequence,
                charge,
                mono_mass_observed,
                intensity_ms1,
                intensity_ms2,
                collision_energy,
                collision_energy_calibrated,
                retention_time,
                retention_time_calibrated,
                retention_time_projected,
                inverse_ion_mobility,
                inverse_ion_mobility_calibrated,
                prosit_predicted_intensities,
                re_score,
                q_value,
                posterior_error_probability,
                py_fragment_intensity_prediction.map(|f| f.inner),
            ),
        }
    }

    #[getter]
    pub fn spec_idx(&self) -> String {
        self.inner.spec_idx.clone()
    }
    
    #[setter]
    pub fn set_spec_idx(&mut self, value: String) {
        self.inner.spec_idx = value;
    }

    #[getter]
    pub fn peptide_idx(&self) -> u32 {
        self.inner.peptide_idx
    }
    
    #[setter]
    pub fn set_peptide_idx(&mut self, value: u32) {
        self.inner.peptide_idx = value;
    }

    #[getter]
    pub fn proteins(&self) -> Vec<String> {
        self.inner.proteins.clone()
    }
    
    #[setter]
    pub fn set_proteins(&mut self, value: Vec<String>) {
        self.inner.proteins = value;
    }

    #[getter]
    pub fn hyperscore(&self) -> f64 {
        self.inner.hyperscore
    }
    
    #[setter]
    pub fn set_hyperscore(&mut self, value: f64) {
        self.inner.hyperscore = value;
    }

    #[getter]
    pub fn sage_feature(&self) -> PyFeature {
        PyFeature {
            inner: self.inner.sage_feature.clone(),
        }
    }
    
    #[setter]
    pub fn set_sage_feature(&mut self, value: PyFeature) {
        self.inner.sage_feature = value.inner.clone();
    }

    #[getter]
    pub fn sequence(&self) -> Option<String> {
        Some(self.inner.clone().sequence?.sequence)
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }
    
    #[setter]
    pub fn set_charge(&mut self, value: Option<u8>) {
        self.inner.charge = value;
    }
    
    #[getter]
    pub fn mono_mass_observed(&self) -> Option<f32> {
        self.inner.mono_mass_observed
    }
    
    #[setter]
    pub fn set_mono_mass_observed(&mut self, value: Option<f32>) {
        self.inner.mono_mass_observed = value;
    }

    #[getter]
    pub fn mono_mass_calculated(&self) -> Option<f32> {
        self.inner.mono_mass_calculated
    }
    
    #[setter]
    pub fn set_mono_mass_calculated(&mut self, value: Option<f32>) {
        self.inner.mono_mass_calculated = value;
    }

    #[getter]
    pub fn mono_mz_calculated(&self) -> Option<f32> {
        self.inner.mono_mz_calculated
    }
    
    #[setter]
    pub fn set_mono_mz_calculated(&mut self, value: Option<f32>) {
        self.inner.mono_mz_calculated = value;
    }
    
    #[getter]
    pub fn intensity_ms1(&self) -> Option<f32> {
        self.inner.intensity_ms1
    }
    
    #[setter]
    pub fn set_intensity_ms1(&mut self, value: Option<f32>) {
        self.inner.intensity_ms1 = value;
    }
    
    #[getter]
    pub fn intensity_ms2(&self) -> Option<f32> {
        self.inner.intensity_ms2
    }
    
    #[setter]
    pub fn set_intensity_ms2(&mut self, value: Option<f32>) {
        self.inner.intensity_ms2 = value;
    }
    
    #[getter]
    pub fn collision_energy(&self) -> Option<f32> {
        self.inner.collision_energy
    }
    
    #[setter]
    pub fn set_collision_energy(&mut self, value: Option<f32>) {
        self.inner.collision_energy = value;
    }
    
    #[getter]
    pub fn collision_energy_calibrated(&self) -> Option<f32> {
        self.inner.collision_energy_calibrated
    }
    
    #[setter]
    pub fn set_collision_energy_calibrated(&mut self, value: Option<f32>) {
        self.inner.collision_energy_calibrated = value;
    }
    
    #[getter]
    pub fn retention_time(&self) -> Option<f32> {
        self.inner.retention_time
    }
    
    #[setter]
    pub fn set_retention_time(&mut self, value: Option<f32>) {
        self.inner.retention_time = value;
    }
    
    #[getter]
    pub fn retention_time_calibrated(&self) -> Option<f32> {
        self.inner.retention_time_calibrated
    }
    
    #[setter]
    pub fn set_retention_time_calibrated(&mut self, value: Option<f32>) {
        self.inner.retention_time_calibrated = value;
    }

    #[getter]
    pub fn retention_time_projected(&self) -> Option<f32> {
        self.inner.retention_time_projected
    }

    #[setter]
    pub fn set_retention_time_projected(&mut self, value: Option<f32>) {
        self.inner.retention_time_projected = value;
    }
    
    #[getter]
    pub fn inverse_ion_mobility(&self) -> Option<f32> {
        self.inner.inverse_ion_mobility
    }
    
    #[setter]
    pub fn set_inverse_ion_mobility(&mut self, value: Option<f32>) {
        self.inner.inverse_ion_mobility = value;
    }
    
    #[getter]
    pub fn inverse_ion_mobility_calibrated(&self) -> Option<f32> {
        self.inner.inverse_ion_mobility_calibrated
    }
    
    #[setter]
    pub fn set_inverse_ion_mobility_calibrated(&mut self, value: Option<f32>) {
        self.inner.inverse_ion_mobility_calibrated = value;
    }
    
    #[getter]
    pub fn prosit_predicted_intensities(&self) -> Option<Vec<f32>> {
        self.inner.prosit_predicted_intensities.clone()
    }
    
    #[setter]
    pub fn set_prosit_predicted_intensities(&mut self, value: Option<Vec<f32>>) {
        self.inner.prosit_predicted_intensities = value;
        self.inner.calculate_fragment_intensity_prediction();
    }
    
    #[getter]
    pub fn re_score(&self) -> Option<f64> {
        self.inner.re_score
    }
    
    #[setter]
    pub fn set_re_score(&mut self, value: Option<f64>) {
        self.inner.re_score = value;
    }
    
    #[getter]
    pub fn q_value(&self) -> Option<f64> {
        self.inner.q_value
    }
    
    #[setter]
    pub fn set_q_value(&mut self, value: Option<f64>) {
        self.inner.q_value = value;
    }
    
    #[getter]
    pub fn posterior_error_probability(&self) -> Option<f64> {
        self.inner.posterior_error_probability
    }
    
    #[setter]
    pub fn set_posterior_error_probability(&mut self, value: Option<f64>) {
        self.inner.posterior_error_probability = value;
    }
    
    #[getter]
    pub fn decoy(&self) -> bool {
        self.inner.decoy
    }
    
    #[setter]
    pub fn set_decoy(&mut self, value: bool) {
        self.inner.decoy = value;
    }
    
    #[getter]
    pub fn external_features(&self) -> BTreeMap<String, f32> {
        self.inner.external_features.clone()
    }
    
    #[setter]
    pub fn set_external_features(&mut self, value: BTreeMap<String, f32>) {
        self.inner.external_features = value;
    }
}

#[pyclass]
#[derive(Clone, Serialize, Deserialize)]
pub struct PyScoreType {
    pub inner: ScoreType,
}

#[pymethods]
impl PyScoreType {
    #[new]
    pub fn new(name: &str) -> Self {
        let score = match name.to_lowercase().as_str() {
            "hyperscore" => SageHyperScore,
            "openmshyperscore" => OpenMSHyperScore,
            _ => panic!("Invalid score type: {}", name),
        };
        
        PyScoreType {
            inner: score
        }
    }

    pub fn to_str(&self) -> String {
        match self.inner {
            SageHyperScore => "hyperscore".to_string(),
            OpenMSHyperScore => "openmshyperscore".to_string(),
        }
    }
}

#[pyclass]
#[derive(Clone, Serialize)]
pub struct PyFragments {
    pub inner: Fragments,
}

#[pymethods]
impl PyFragments {
    #[new]
    pub fn new(
        charges: Vec<i32>,
        kinds: Vec<PyKind>,
        fragment_ordinals: Vec<i32>,
        intensities: Vec<f32>,
        mz_calculated: Vec<f32>,
        mz_experimental: Vec<f32>,
    ) -> Self {
        PyFragments {
            inner: Fragments {
                charges,
                kinds: kinds.into_iter().map(|k| k.inner).collect(),
                fragment_ordinals,
                intensities,
                mz_calculated,
                mz_experimental,
            },
        }
    }

    #[getter]
    pub fn charges(&self) -> Vec<i32> {
        self.inner.charges.clone()
    }

    #[getter]
    pub fn kinds(&self) -> Vec<PyKind> {
        self.inner.kinds.iter().map(|k| PyKind { inner: *k }).collect()
    }

    #[getter]
    pub fn fragment_ordinals(&self) -> Vec<i32> {
        self.inner.fragment_ordinals.clone()
    }

    #[getter]
    pub fn intensities(&self) -> Vec<f32> {
        self.inner.intensities.clone()
    }

    #[getter]
    pub fn mz_calculated(&self) -> Vec<f32> {
        self.inner.mz_calculated.clone()
    }

    #[getter]
    pub fn mz_experimental(&self) -> Vec<f32> {
        self.inner.mz_experimental.clone()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyFeature {
    pub inner: Feature,
}

#[pymethods]
impl PyFeature {
    #[new]
    pub fn new(
        peptide_idx: PyPeptideIx,
        psm_id: usize,
        peptide_len: usize,
        spec_id: String,
        file_id: usize,
        rank: u32,
        label: i32,
        expmass: f32,
        calcmass: f32,
        charge: u8,
        delta_mass: f32,
        isotope_error: f32,
        average_ppm: f32,
        hyperscore: f64,
        delta_next: f64,
        delta_best: f64,
        matched_peaks: u32,
        longest_b: u32,
        longest_y: u32,
        longest_y_pct: f32,
        missed_cleavages: u8,
        matched_intensity_pct: f32,
        scored_candidates: u32,
        poisson: f64,
        discriminant_score: f32,
        posterior_error: f32,
        spectrum_q: f32,
        peptide_q: f32,
        protein_q: f32,
        ms2_intensity: f32,
        rt: Option<f32>,
        aligned_rt: Option<f32>,
        predicted_rt: Option<f32>,
        delta_rt_model: Option<f32>,
        ims: Option<f32>,
        predicted_ims: Option<f32>,
        delta_ims_model: Option<f32>,
        fragments: Option<PyFragments>,
    ) -> Self {
        PyFeature {
            inner: Feature {
                peptide_idx: peptide_idx.inner,
                psm_id,
                peptide_len,
                spec_id,
                file_id,
                rank,
                label,
                expmass,
                calcmass,
                charge,
                rt: rt.unwrap_or(0.0),
                aligned_rt: aligned_rt.unwrap_or(0.0),
                predicted_rt: predicted_rt.unwrap_or(0.0),
                delta_rt_model: delta_rt_model.unwrap_or(0.0),
                ims: ims.unwrap_or(0.0),
                predicted_ims: predicted_ims.unwrap_or(0.0),
                delta_ims_model: delta_ims_model.unwrap_or(0.0),
                delta_mass,
                isotope_error,
                average_ppm,
                hyperscore,
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
                discriminant_score,
                posterior_error,
                spectrum_q,
                peptide_q,
                protein_q,
                ms2_intensity,
                fragments: fragments.map(|f| f.inner),
            },
        }
    }

    #[getter]
    pub fn peptide_idx(&self) -> PyPeptideIx {
        PyPeptideIx {
            inner: self.inner.peptide_idx,
        }
    }

    #[getter]
    pub fn psm_id(&self) -> usize {
        self.inner.psm_id
    }

    #[getter]
    pub fn peptide_len(&self) -> usize {
        self.inner.peptide_len
    }

    #[getter]
    pub fn spec_id(&self) -> String {
        self.inner.spec_id.clone()
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn rank(&self) -> u32 {
        self.inner.rank
    }

    #[getter]
    pub fn label(&self) -> i32 {
        self.inner.label
    }

    #[getter]
    pub fn expmass(&self) -> f32 {
        self.inner.expmass
    }

    #[getter]
    pub fn calcmass(&self) -> f32 {
        self.inner.calcmass
    }

    #[getter]
    pub fn charge(&self) -> u8 {
        self.inner.charge
    }

    #[getter]
    pub fn rt(&self) -> f32 {
        self.inner.rt
    }

    #[getter]
    pub fn aligned_rt(&self) -> f32 {
        self.inner.aligned_rt
    }

    #[getter]
    pub fn predicted_rt(&self) -> f32 {
        self.inner.predicted_rt
    }

    #[getter]
    pub fn delta_rt_model(&self) -> f32 {
        self.inner.delta_rt_model
    }

    #[getter]
    pub fn delta_mass(&self) -> f32 {
        self.inner.delta_mass
    }

    #[getter]
    pub fn isotope_error(&self) -> f32 {
        self.inner.isotope_error
    }

    #[getter]
    pub fn average_ppm(&self) -> f32 {
        self.inner.average_ppm
    }

    #[getter]
    pub fn hyperscore(&self) -> f64 {
        self.inner.hyperscore
    }

    #[getter]
    pub fn delta_next(&self) -> f64 {
        self.inner.delta_next
    }

    #[getter]
    pub fn delta_best(&self) -> f64 {
        self.inner.delta_best
    }

    #[getter]
    pub fn matched_peaks(&self) -> u32 {
        self.inner.matched_peaks
    }

    #[getter]
    pub fn longest_b(&self) -> u32 {
        self.inner.longest_b
    }

    #[getter]
    pub fn longest_y(&self) -> u32 {
        self.inner.longest_y
    }

    #[getter]
    pub fn longest_y_pct(&self) -> f32 {
        self.inner.longest_y_pct
    }

    #[getter]
    pub fn missed_cleavages(&self) -> u8 {
        self.inner.missed_cleavages
    }

    #[getter]
    pub fn matched_intensity_pct(&self) -> f32 {
        self.inner.matched_intensity_pct
    }

    #[getter]
    pub fn scored_candidates(&self) -> u32 {
        self.inner.scored_candidates
    }

    #[getter]
    pub fn poisson(&self) -> f64 {
        self.inner.poisson
    }

    #[getter]
    pub fn discriminant_score(&self) -> f32 {
        self.inner.discriminant_score
    }

    #[getter]
    pub fn posterior_error(&self) -> f32 {
        self.inner.posterior_error
    }

    #[getter]
    pub fn spectrum_q(&self) -> f32 {
        self.inner.spectrum_q
    }

    #[getter]
    pub fn peptide_q(&self) -> f32 {
        self.inner.peptide_q
    }

    #[getter]
    pub fn protein_q(&self) -> f32 {
        self.inner.protein_q
    }

    #[getter]
    pub fn ms2_intensity(&self) -> f32 {
        self.inner.ms2_intensity
    }

    #[getter]
    pub fn fragments(&self) -> Option<PyFragments> {
        self.inner.fragments.as_ref().map(|f| PyFragments { inner: f.clone() })
    }

    #[getter]
    pub fn ims(&self) -> f32 {
        self.inner.ims
    }

    #[getter]
    pub fn predicted_ims(&self) -> f32 {
        self.inner.predicted_ims
    }

    #[getter]
    pub fn delta_ims_model(&self) -> f32 {
        self.inner.delta_ims_model
    }
}

#[pyclass]
pub struct PyScorer {
    pub precursor_tolerance: PyTolerance,
    pub fragment_tolerance: PyTolerance,
    pub min_matched_peaks: u16,
    pub min_isotope_err: i8,
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub report_psms: usize,
    pub wide_window: bool,
    pub annotate_matches: bool,
    pub override_precursor_charge: bool,
    pub expected_mods: HashSet<String>,
    pub score_type: Option<PyScoreType>,
}

#[pymethods]
impl PyScorer {
    #[new]
    pub fn new(
        precursor_tolerance: PyTolerance,
        fragment_tolerance: PyTolerance,
        min_matched_peaks: u16,
        min_isotope_err: i8,
        max_isotope_err: i8,
        min_precursor_charge: u8,
        max_precursor_charge: u8,
        chimera: bool,
        report_psms: usize,
        wide_window: bool,
        annotate_matches: bool,
        override_precursor_charge: bool,
        expected_mods: HashSet<String>,
        max_fragment_charge: Option<u8>,
        score_type: Option<PyScoreType>,
    ) -> Self {
        PyScorer {
            precursor_tolerance,
            fragment_tolerance,
            min_matched_peaks,
            min_isotope_err,
            max_isotope_err,
            min_precursor_charge,
            max_precursor_charge,
            max_fragment_charge,
            chimera,
            report_psms,
            wide_window,
            annotate_matches,
            override_precursor_charge,
            score_type,
            expected_mods,
        }
    }

    pub fn score(&self, db: &PyIndexedDatabase, spectrum: &PyProcessedSpectrum) -> Vec<PyFeature> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().unwrap().inner,
        };
        let features = scorer.score(&spectrum.inner);
        features
            .into_iter()
            .map(|f| PyFeature { inner: f })
            .collect()
    }

    pub fn score_collection(
        &self,
        db: &PyIndexedDatabase,
        spectra: Vec<PyProcessedSpectrum>,
        num_threads: usize,
    ) -> Vec<Vec<PyFeature>> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().unwrap().inner,
        };
        // Configure the global thread pool to the desired number of threads
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let result = pool.install(|| {
            spectra
                .par_iter()
                .map(|spectrum| {
                    let features = scorer.score(&spectrum.inner);
                    features
                        .into_iter()
                        .map(|f| PyFeature { inner: f })
                        .collect()
                })
                .collect()
        });

        result
    }
    
    pub fn score_candidates(
        &self,
        db: &PyIndexedDatabase,
        spectra: Vec<PyProcessedSpectrum>,
        num_threads: usize,
    ) -> BTreeMap<String, Vec<PyPsm>> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().unwrap().inner,
        };

        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let result: Vec<Vec<Feature>> = pool.install(|| {
            spectra
                .par_iter()
                .map(|spectrum| scorer.score(&spectrum.inner))
                .collect()
        });
        
        let psm_map: BTreeMap<String, Vec<Psm>> = pool.install(|| {
            
            spectra.par_iter().zip(result.into_par_iter())
                
                .map(|(spectrum, features)| {
                    
                    let mut psms = Vec::new();
                    
                    for feature in &features {
                        
                        let peptide = &db.inner[feature.peptide_idx];
                        let decoy = peptide.decoy;
                        let score = feature.hyperscore;
                        
                        let intensity_ms1: f32 = spectrum.inner.precursors.iter().map(|p| p.intensity.unwrap_or(0.0)).sum();
                        let intensity_ms2: f32 = feature.ms2_intensity;
                        let charge = feature.charge;
                        let proteins: Vec<String> = peptide.proteins.iter().map(|arc| (**arc).clone()).collect();
                        
                        let sequence = std::str::from_utf8(&peptide.sequence).unwrap().to_string();
                        let peptide_sequence = sage_sequence_to_unimod_sequence(sequence, &peptide.modifications, &self.expected_mods);
                        
                        let collision_energy = spectrum.collision_energies.first().unwrap_or(&None).unwrap_or(0.0f32);
                        
                        let psm = Psm::new(
                            spectrum.inner.id.clone(),
                            feature.clone().peptide_idx.0,
                            proteins,
                            score,
                            decoy,
                            feature.clone(),
                            Some(peptide_sequence),
                            Some(charge),
                            Some(feature.expmass),
                            Some(intensity_ms1),
                            Some(intensity_ms2),
                            Some(collision_energy),
                            None, // collision_energy_calibrated
                            Some(feature.rt),
                            None, // retention_time_calibrated
                            None, // rt projected
                            Some(feature.ims),
                            None, // inverse_ion_mobility_calibrated
                            None, // prosit_predicted_intensities
                            None, // re_score
                            None, // q_value
                            None, // posterior_error_probability
                            None, // fragment_intensity_prediction
                        );
                        psms.push(psm);
                    }
                    (spectrum.inner.id.clone(), psms)
                })
                .collect()
        });
        
        let mut result: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();
        
        for (spec_id, psms) in psm_map {
            result.insert(spec_id, psms.into_iter().map(|psm| PyPsm { inner: psm }).collect());
        }
        
        // sort by spec_id, hyper_score, peptide_idx, decoy
        for (_, psms) in result.iter_mut() {
            psms.sort_by(|a, b| {
                let a = &a.inner;
                let b = &b.inner;
                a.hyperscore.partial_cmp(&b.hyperscore).unwrap()
                    .then_with(|| a.peptide_idx.partial_cmp(&b.peptide_idx).unwrap())
                    .then_with(|| a.decoy.partial_cmp(&b.decoy).unwrap())
            });
        }
        
        result
    }

    pub fn score_collection_to_psm_collection(
        &self,
        db: &PyIndexedDatabase,
        spectra: Vec<PyProcessedSpectrum>,
        num_threads: usize,
    ) -> BTreeMap<String, Vec<PyPeptideSpectrumMatch>> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().unwrap().inner,
        };

        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let result: Vec<Vec<Feature>> = pool.install(|| {
            spectra
                .par_iter()
                .map(|spectrum| scorer.score(&spectrum.inner))
                .collect()
        });

        let psm_map: BTreeMap<String, Vec<(PeptideSpectrumMatch, Option<Fragments>)>> = pool.install(|| {
            spectra.par_iter().zip(result.into_par_iter())
                .map(|(spectrum, features)| {
                    let mut psms = Vec::new();

                    let mut hash_set: HashSet<(u32, bool)> = HashSet::new();

                    for feature in features {
                        let peptide = &db.inner[feature.peptide_idx];
                        let decoy = peptide.decoy;
                        let score = feature.hyperscore;
                        let intensity_ms1: f32 = spectrum.inner.precursors.iter().map(|p| p.intensity.unwrap_or(0.0)).sum();
                        let intensity_ms2: f32 = feature.ms2_intensity;
                        let charge = feature.charge;
                        let proteins: Vec<String> = peptide.proteins.iter().map(|arc| (**arc).clone()).collect();
                        let sequence = std::str::from_utf8(&peptide.sequence).unwrap().to_string();
                        let fragments = feature.fragments;
                        let collision_energy = spectrum.collision_energies.first().unwrap_or(&None).unwrap_or(0.0f32);

                        let key = (feature.peptide_idx.0, decoy);

                        if hash_set.contains(&key) {
                            continue;
                        }

                        hash_set.insert(key);
                        let maybe_charges = fragments.clone().map(|f| f.charges);
                        let maybe_kinds = fragments.clone().map(|f| f.kinds.iter().map(|k| kind_to_string(*k)).collect());
                        let maybe_fragment_ordinals = fragments.clone().map(|f| f.fragment_ordinals);
                        let maybe_intensities = fragments.clone().map(|f| f.intensities);
                        let maybe_mz_calculated = fragments.clone().map(|f| f.mz_calculated);
                        let maybe_mz_experimental = fragments.clone().map(|f| f.mz_experimental);

                        let psm = PeptideSpectrumMatch::new(
                            spectrum.inner.id.clone(),
                            feature.peptide_idx.0,
                            proteins,
                            decoy,
                            score,
                            feature.rank,
                            Some(feature.expmass),
                            Some(feature.isotope_error),
                            Some(feature.average_ppm),
                            Some(feature.delta_next),
                            Some(feature.delta_best),
                            Some(feature.matched_peaks),
                            Some(feature.longest_b),
                            Some(feature.longest_y),
                            Some(feature.longest_y_pct),
                            Some(feature.missed_cleavages),
                            Some(feature.matched_intensity_pct),
                            Some(feature.scored_candidates),
                            Some(feature.poisson),
                            Some(sage_sequence_to_unimod_sequence(sequence, &peptide.modifications, &self.expected_mods)),
                            Some(charge),
                            Some(feature.rt),
                            None,
                            Some(feature.ims),
                            None,
                            Some(intensity_ms1),
                            Some(intensity_ms2),
                            None,
                            Some(collision_energy as f64),
                            None,
                            None,
                            None,
                            None,
                            maybe_charges,
                            maybe_kinds,
                            maybe_fragment_ordinals,
                            maybe_intensities,
                            maybe_mz_calculated,
                            maybe_mz_experimental,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                            None,
                        );
                        psms.push((psm, fragments));
                    }
                    (spectrum.inner.id.clone(), psms)
                })
                .collect()
        });

        let mut result: Vec<PyPeptideSpectrumMatch> = Vec::new();

        for (_, psms) in psm_map {
            for (psm, fragments) in psms {
                result.push(PyPeptideSpectrumMatch {
                    inner: psm,
                    fragments_observed: fragments.map(|f| PyFragments { inner: f }),
                    fragments_predicted: None,
                });
            }
        }

        // sort by spec_id, hyper_score, peptide_idx, decoy
        result.sort_by(|a, b| {
            let a = &a.inner;
            let b = &b.inner;
            a.spec_idx.cmp(&b.spec_idx)
                .then(b.hyper_score.partial_cmp(&a.hyper_score).unwrap())
                .then(a.peptide_idx.cmp(&b.peptide_idx))
                .then(a.decoy.cmp(&b.decoy))
        });

        let mut ret_map: BTreeMap<String, Vec<PyPeptideSpectrumMatch>> = BTreeMap::new();

        for psm in result {
            let key = psm.inner.spec_idx.clone();
            if !ret_map.contains_key(&key) {
                ret_map.insert(key.clone(), Vec::new());
            }
            ret_map.get_mut(&key).unwrap().push(psm);
        }

        ret_map
    }

    pub fn score_chimera_fast(
        &self,
        db: &PyIndexedDatabase,
        query: &PyProcessedSpectrum,
    ) -> Vec<PyFeature> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().unwrap().inner,
        };
        let features = scorer.score_chimera_fast(&query.inner);
        features
            .into_iter()
            .map(|f| PyFeature { inner: f })
            .collect()
    }

    pub fn score_standard(
        &self,
        db: &PyIndexedDatabase,
        query: &PyProcessedSpectrum,
    ) -> Vec<PyFeature> {
        let scorer = Scorer {
            db: &db.inner,
            precursor_tol: self.precursor_tolerance.inner.clone(),
            fragment_tol: self.fragment_tolerance.inner.clone(),
            min_matched_peaks: self.min_matched_peaks,
            min_isotope_err: self.min_isotope_err,
            max_isotope_err: self.max_isotope_err,
            min_precursor_charge: self.min_precursor_charge,
            max_precursor_charge: self.max_precursor_charge,
            max_fragment_charge: self.max_fragment_charge,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
            override_precursor_charge: self.override_precursor_charge,
            score_type: self.score_type.clone().map(|s| s.inner).unwrap(),
        };
        let features = scorer.score_standard(&query.inner);
        features
            .into_iter()
            .map(|f| PyFeature { inner: f })
            .collect()
    }

    #[getter]
    pub fn precursor_tolerance(&self) -> PyTolerance {
        self.precursor_tolerance.clone()
    }

    #[getter]
    pub fn fragment_tolerance(&self) -> PyTolerance {
        self.fragment_tolerance.clone()
    }

    #[getter]
    pub fn min_matched_peaks(&self) -> u16 {
        self.min_matched_peaks
    }

    #[getter]
    pub fn min_isotope_err(&self) -> i8 {
        self.min_isotope_err
    }

    #[getter]
    pub fn max_isotope_err(&self) -> i8 {
        self.max_isotope_err
    }

    #[getter]
    pub fn min_precursor_charge(&self) -> u8 {
        self.min_precursor_charge
    }

    #[getter]
    pub fn max_precursor_charge(&self) -> u8 {
        self.max_precursor_charge
    }

    #[getter]
    pub fn max_fragment_charge(&self) -> Option<u8> {
        self.max_fragment_charge
    }

    #[getter]
    pub fn chimera(&self) -> bool {
        self.chimera
    }

    #[getter]
    pub fn report_psms(&self) -> usize {
        self.report_psms
    }

    #[getter]
    pub fn wide_window(&self) -> bool {
        self.wide_window
    }
}

#[pyclass]
#[derive(Clone, Serialize)]
pub struct PyPeptideSpectrumMatch {
    pub inner: PeptideSpectrumMatch,
    pub fragments_observed: Option<PyFragments>,
    pub fragments_predicted: Option<PyFragments>,
}

#[pymethods]
impl PyPeptideSpectrumMatch {
    #[new]
    fn new(
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
        fragments_observed: Option<PyFragments>,
        fragments_predicted: Option<PyFragments>,
        re_score: Option<f64>,
        cosine_similarity: Option<f32>,
        file_name: Option<String>,
        mz_calibration_ppm: Option<f32>,
        projected_rt: Option<f32>,
        beta_score: Option<f64>,
        posterior_error_prob: Option<f64>,
        prosit_intensities: Option<Vec<f32>>,
        spectral_entropy_similarity: Option<f32>,
        spectral_correlation_similarity_pearson: Option<f32>,
        spectral_correlation_similarity_spearman: Option<f32>,
        spectral_normalized_intensity_difference: Option<f32>,
    ) -> Self {

        let maybe_charges = fragments_observed.clone().map(|f| f.inner.charges);
        let maybe_kinds = fragments_observed.clone().map(|f| f.inner.kinds.iter().map(|k| kind_to_string(*k)).collect());
        let maybe_fragment_ordinals = fragments_observed.clone().map(|f| f.inner.fragment_ordinals);
        let maybe_intensities = fragments_observed.clone().map(|f| f.inner.intensities);
        let maybe_mz_calculated = fragments_observed.clone().map(|f| f.inner.mz_calculated);
        let maybe_mz_experimental = fragments_observed.clone().map(|f| f.inner.mz_experimental);

        PyPeptideSpectrumMatch {
            inner: PeptideSpectrumMatch::new(
                spec_idx,
                peptide_idx,
                proteins,
                decoy,
                hyper_score,
                rank,
                mono_mass_observed,
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
                sequence,
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
                maybe_charges,
                maybe_kinds,
                maybe_fragment_ordinals,
                maybe_intensities,
                maybe_mz_calculated,
                maybe_mz_experimental,
                mz_calibration_ppm,
                projected_rt,
                beta_score,
                posterior_error_prob,
                prosit_intensities,
                spectral_entropy_similarity,
                spectral_correlation_similarity_pearson,
                spectral_correlation_similarity_spearman,
                spectral_normalized_intensity_difference,
            ),
            fragments_observed,
            fragments_predicted,
        }
    }

    #[getter]
    pub fn spec_idx(&self) -> &str {
        &self.inner.spec_idx
    }

    #[getter]
    pub fn peptide_idx(&self) -> u32 {
        self.inner.peptide_idx
    }

    #[getter]
    pub fn proteins(&self) -> Vec<String> {
        self.inner.proteins.clone()
    }

    #[getter]
    pub fn decoy(&self) -> bool {
        self.inner.decoy
    }

    #[getter]
    pub fn hyper_score(&self) -> f64 {
        self.inner.hyper_score
    }

    #[setter]
    pub fn set_hyper_score(&mut self, hyper_score: f64) {
        self.inner.hyper_score = hyper_score;
    }

    #[getter]
    pub fn re_score(&self) -> Option<f64> {
        self.inner.re_score
    }

    #[setter]
    pub fn set_re_score(&mut self, re_score: f64) {
        self.inner.re_score = Some(re_score);
    }

    #[getter]
    pub fn rank(&self) -> u32 {
        self.inner.rank
    }

    #[getter]
    pub fn mono_mz_calculated(&self) -> Option<f32> {
        self.inner.mono_mz_calculated
    }

    #[getter]
    pub fn mono_mass_calculated(&self) -> Option<f32> {
        self.inner.mono_mass_calculated
    }

    #[getter]
    pub fn mono_mass_observed(&self) -> Option<f32> {
        self.inner.mono_mass_observed
    }

    #[getter]
    pub fn isotope_error(&self) -> Option<f32> {
        self.inner.isotope_error
    }

    #[getter]
    pub fn average_ppm(&self) -> Option<f32> {
        self.inner.average_ppm
    }

    #[getter]
    pub fn delta_next(&self) -> Option<f64> {
        self.inner.delta_next
    }

    #[getter]
    pub fn delta_best(&self) -> Option<f64> {
        self.inner.delta_best
    }

    #[getter]
    pub fn matched_peaks(&self) -> Option<u32> {
        self.inner.matched_peaks
    }

    #[getter]
    pub fn longest_b(&self) -> Option<u32> {
        self.inner.longest_b
    }

    #[getter]
    pub fn longest_y(&self) -> Option<u32> {
        self.inner.longest_y
    }

    #[getter]
    pub fn longest_y_pct(&self) -> Option<f32> {
        self.inner.longest_y_pct
    }

    #[getter]
    pub fn missed_cleavages(&self) -> Option<u8> {
        self.inner.missed_cleavages
    }

    #[getter]
    pub fn matched_intensity_pct(&self) -> Option<f32> {
        self.inner.matched_intensity_pct
    }

    #[getter]
    pub fn scored_candidates(&self) -> Option<u32> {
        self.inner.scored_candidates
    }

    #[getter]
    pub fn poisson(&self) -> Option<f64> {
        self.inner.poisson
    }

    #[getter]
    pub fn peptide_sequence(&self) -> Option<String> {
        match self.inner.peptide_sequence {
            Some(ref seq) => Some(seq.sequence.clone()),
            None => None,
        }
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }
    #[getter]
    pub fn fragments_observed(&mut self) -> Option<PyFragments> {
        if self.fragments_observed.is_none() {
            if self.inner.fragment_charges.is_none() {
                return None;
            } else {
                let charges = self.inner.fragment_charges.as_ref().unwrap();
                let kinds = self.inner.fragment_ion_types.as_ref().unwrap();
                let fragment_ordinals = self.inner.fragment_ordinals.as_ref().unwrap();
                let intensities = self.inner.fragment_intensities.as_ref().unwrap();
                let mz_calculated = self.inner.fragment_mz_calculated.as_ref().unwrap();
                let mz_experimental = self.inner.fragment_mz_observed.as_ref().unwrap();

                let fragments = Fragments {
                    charges: charges.clone(),
                    kinds: kinds.into_iter().map(|k| string_to_kind(k)).collect(),
                    fragment_ordinals: fragment_ordinals.clone(),
                    intensities: intensities.clone(),
                    mz_calculated: mz_calculated.clone(),
                    mz_experimental: mz_experimental.clone(),
                };

                self.fragments_observed = Some(PyFragments { inner: fragments });
            }
        }

        self.fragments_observed.clone()
    }

    #[setter]
    pub fn set_fragments_observed(&mut self, fragments: Option<PyFragments>) {
        self.fragments_observed = fragments;
    }

    #[getter]
    pub fn fragments_predicted(&self) -> Option<PyFragments> {
        self.fragments_predicted.clone()
    }

    #[setter]
    pub fn set_fragments_predicted(&mut self, fragments: Option<PyFragments>) {
        self.fragments_predicted = fragments;
    }

    #[getter]
    pub fn retention_time_observed(&self) -> Option<f32> {
        self.inner.retention_time_observed
    }

    #[getter]
    pub fn retention_time_predicted(&self) -> Option<f32> {
        self.inner.retention_time_predicted
    }

    #[setter]
    pub fn set_retention_time_predicted(&mut self, retention_time_predicted: f32) {
        self.inner.retention_time_predicted = Some(retention_time_predicted);
    }

    #[getter]
    pub fn inverse_mobility_observed(&self) -> Option<f32> {
        self.inner.inverse_mobility_observed
    }

    #[getter]
    pub fn inverse_mobility_predicted(&self) -> Option<f32> {
        self.inner.inverse_mobility_predicted
    }

    #[setter]
    pub fn set_inverse_mobility_predicted(&mut self, inverse_mobility_predicted: f32) {
        self.inner.inverse_mobility_predicted = Some(inverse_mobility_predicted);
    }

    #[getter]
    pub fn intensity_ms1(&self) -> Option<f32> {
        self.inner.intensity_ms1
    }

    #[getter]
    pub fn intensity_ms2(&self) -> Option<f32> {
        self.inner.intensity_ms2
    }

    #[getter]
    pub fn q_value(&self) -> Option<f64> {
        self.inner.q_value
    }

    #[setter]
    pub fn set_q_value(&mut self, q_value: f64) {
        self.inner.q_value = Some(q_value);
    }

    #[getter]
    pub fn collision_energy(&self) -> Option<f64> {
        self.inner.collision_energy
    }

    #[setter]
    pub fn set_collision_energy(&mut self, collision_energy: f64) {
        self.inner.collision_energy = Some(collision_energy);
    }

    #[getter]
    pub fn collision_energy_calibrated(&self) -> Option<f64> {
        self.inner.collision_energy_calibrated
    }

    #[setter]
    pub fn set_collision_energy_calibrated(&mut self, collision_energy_calibrated: f64) {
        self.inner.collision_energy_calibrated = Some(collision_energy_calibrated);
    }

    #[getter]
    pub fn cosine_similarity(&self) -> Option<f32> {
        self.inner.cosine_similarity
    }

    #[setter]
    pub fn set_cosine_similarity(&mut self, cosine_similarity: Option<f32>) {
        self.inner.cosine_similarity = cosine_similarity;
    }

    #[getter]
    pub fn file_name(&self) -> Option<String> {
        self.inner.file_name.clone()
    }

    #[setter]
    pub fn set_file_name(&mut self, file_name: String) {
        self.inner.file_name = Some(file_name);
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(&self.inner).unwrap()
    }

    #[getter]
    pub fn mz_calibration_ppm(&self) -> Option<f32> {
        self.inner.mz_calibration_ppm
    }

    #[setter]
    pub fn set_mz_calibration_ppm(&mut self, mz_calibration_ppm: f32) {
        self.inner.mz_calibration_ppm = Some(mz_calibration_ppm);
    }

    #[getter]
    pub fn projected_rt(&self) -> Option<f32> {
        self.inner.projected_rt
    }

    #[setter]
    pub fn set_projected_rt(&mut self, projected_rt: f32) {
        self.inner.projected_rt = Some(projected_rt);
    }

    #[getter]
    pub fn beta_score(&self) -> Option<f64> {
        self.inner.beta_score
    }

    #[setter]
    pub fn set_beta_score(&mut self, beta_score: f64) {
        self.inner.beta_score = Some(beta_score);
    }

    #[getter]
    pub fn posterior_error_prob(&self) -> Option<f64> {
        self.inner.posterior_error_prob
    }

    #[setter]
    pub fn set_posterior_error_prob(&mut self, posterior_error_prob: f64) {
        self.inner.posterior_error_prob = Some(posterior_error_prob);
    }

    #[getter]
    pub fn prosit_intensities(&self) -> Option<Vec<f32>> {
        self.inner.prosit_intensities.clone()
    }

    #[setter]
    pub fn set_prosit_intensities(&mut self, prosit_intensities: Vec<f32>) {
        self.inner.prosit_intensities = Some(prosit_intensities);
    }

    #[getter]
    pub fn spectral_entropy_similarity(&self) -> Option<f32> {
        self.inner.spectral_entropy_similarity
    }

    #[setter]
    pub fn set_spectral_entropy_similarity(&mut self, spectral_entropy_similarity: f32) {
        self.inner.spectral_entropy_similarity = Some(spectral_entropy_similarity);
    }

    #[getter]
    pub fn spectral_correlation_similarity_pearson(&self) -> Option<f32> {
        self.inner.spectral_correlation_similarity_pearson
    }

    #[setter]
    pub fn set_spectral_correlation_similarity_pearson(&mut self, spectral_correlation_similarity_pearson: f32) {
        self.inner.spectral_correlation_similarity_pearson = Some(spectral_correlation_similarity_pearson);
    }

    #[getter]
    pub fn spectral_correlation_similarity_spearman(&self) -> Option<f32> {
        self.inner.spectral_correlation_similarity_spearman
    }

    #[setter]
    pub fn set_spectral_correlation_similarity_spearman(&mut self, spectral_correlation_similarity_spearman: f32) {
        self.inner.spectral_correlation_similarity_spearman = Some(spectral_correlation_similarity_spearman);
    }

    #[getter]
    pub fn spectral_normalized_intensity_difference(&self) -> Option<f32> {
        self.inner.spectral_normalized_intensity_difference
    }

    #[setter]
    pub fn set_spectral_normalized_intensity_difference(&mut self, spectral_normalized_intensity_difference: f32) {
        self.inner.spectral_normalized_intensity_difference = Some(spectral_normalized_intensity_difference);
    }

    pub fn prosit_fragment_map(&self) -> Option<BTreeMap<(u32, i32, i32), f32>> {
        self.inner.prosit_intensities.clone().map(|intensities| {
            flat_prosit_array_to_fragments_map(intensities)
        })
    }

    pub fn observed_fragment_map(&self) -> Option<BTreeMap<(u32, i32, i32), f32>> {
        self.fragments_observed.clone().map(|f| py_fragments_to_fragments_map(&f, true))
    }
}

#[pyfunction]
pub fn psm_from_json(json: &str) -> PyPeptideSpectrumMatch {
    let psm: PeptideSpectrumMatch = serde_json::from_str(json).unwrap();
    PyPeptideSpectrumMatch {
        inner: psm,
        fragments_observed: None,
        fragments_predicted: None,
    }
}

#[pyfunction]
pub fn peptide_spectrum_match_from_json(json: &str) -> PyPsm {
    let psm: Psm = serde_json::from_str(json).unwrap();
    PyPsm {
        inner: psm
    }
}

#[pyfunction]
pub fn prosit_intensities_to_py_fragments(
    flat_intensities: Vec<f32>,
) -> PyFragments {
    let fragments = flat_prosit_array_to_fragments_map(flat_intensities);

    let mut predicted_kinds_b: Vec<Kind> = Vec::new();
    let mut predicted_kinds_y: Vec<Kind> = Vec::new();

    let mut predicted_fragment_ordinals_b = Vec::new();
    let mut predicted_fragment_ordinals_y = Vec::new();

    let mut predicted_charges_b = Vec::new();
    let mut predicted_charges_y = Vec::new();

    let mut predicted_intensities_b = Vec::new();
    let mut predicted_intensities_y = Vec::new();

    for (key, _) in fragments.iter() {

        let (kind, charge, fragment_ordinal) = key;

        let predicted_intensity = fragments.get(key).unwrap_or(&0.0);

        let kind = match kind {
            0 => Kind::B,
            1 => Kind::Y,
            _ => panic!("Invalid kind"),
        };

        if kind == Kind::B {
            predicted_kinds_b.push(kind);
            predicted_charges_b.push(*charge);
            predicted_fragment_ordinals_b.push(*fragment_ordinal);
            predicted_intensities_b.push(*predicted_intensity);
        } else {
            predicted_kinds_y.push(kind);
            predicted_charges_y.push(*charge);
            predicted_fragment_ordinals_y.push(*fragment_ordinal);
            predicted_intensities_y.push(*predicted_intensity);
        }
    }

    // invert the order of y fragments
    predicted_kinds_y.reverse();
    predicted_charges_y.reverse();
    predicted_fragment_ordinals_y.reverse();
    predicted_intensities_y.reverse();

    let fragments = Fragments {
        charges: predicted_charges_b.iter().chain(predicted_charges_y.iter()).cloned().collect(),
        kinds: predicted_kinds_b.iter().chain(predicted_kinds_y.iter()).cloned().collect(),
        fragment_ordinals: predicted_fragment_ordinals_b.iter().chain(predicted_fragment_ordinals_y.iter()).cloned().collect(),
        intensities: predicted_intensities_b.iter().chain(predicted_intensities_y.iter()).cloned().collect(),
        mz_calculated: Vec::new(),
        mz_experimental: Vec::new(),
    };

    PyFragments {
        inner: fragments
    }
}

#[pyfunction]
pub fn prosit_intensities_to_py_fragments_par(
    flat_intensities: Vec<Vec<f32>>,
    num_threads: usize
) -> Vec<PyFragments> {
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let result = pool.install(|| {
        flat_intensities.par_iter()
            .map(|intensities| {
                prosit_intensities_to_py_fragments(intensities.clone())
            })
            .collect()
    });

    result
}

#[pyfunction]
pub fn associate_psm_with_prosit_predicted_intensities(
    psm: PyPeptideSpectrumMatch,
    flat_intensities: Vec<f32>,
) -> PyPeptideSpectrumMatch {

    let intensity_copy = flat_intensities.clone();

    let fragments_observed = &psm.fragments_observed.unwrap();
    let observed_map = py_fragments_to_fragments_map(fragments_observed, true);
    let predicted_map = flat_prosit_array_to_fragments_map(flat_intensities);

    let mut predicted_kinds_b: Vec<Kind> = Vec::new();
    let mut predicted_kinds_y: Vec<Kind> = Vec::new();

    let mut predicted_fragment_ordinals_b = Vec::new();
    let mut predicted_fragment_ordinals_y = Vec::new();

    let mut predicted_charges_b = Vec::new();
    let mut predicted_charges_y = Vec::new();

    let mut predicted_intensities_b = Vec::new();
    let mut predicted_intensities_y = Vec::new();

    for (key, _) in observed_map.iter() {

        let (kind, charge, fragment_ordinal) = key;

        let predicted_intensity = predicted_map.get(key).unwrap_or(&0.0);

        let kind = match kind {
            0 => Kind::B,
            1 => Kind::Y,
            _ => panic!("Invalid kind"),
        };

        if kind == Kind::B {
            predicted_kinds_b.push(kind);
            predicted_charges_b.push(*charge);
            predicted_fragment_ordinals_b.push(*fragment_ordinal);
            predicted_intensities_b.push(*predicted_intensity);
        } else {
            predicted_kinds_y.push(kind);
            predicted_charges_y.push(*charge);
            predicted_fragment_ordinals_y.push(*fragment_ordinal);
            predicted_intensities_y.push(*predicted_intensity);
        }
    }

    // invert the order of y fragments
    predicted_kinds_y.reverse();
    predicted_charges_y.reverse();
    predicted_fragment_ordinals_y.reverse();
    predicted_intensities_y.reverse();

    let fragments_predicted = PyFragments {
        inner: Fragments {
            // concat the two lists
            charges: predicted_charges_b.iter().chain(predicted_charges_y.iter()).cloned().collect(),
            kinds: predicted_kinds_b.iter().chain(predicted_kinds_y.iter()).cloned().collect(),
            fragment_ordinals: predicted_fragment_ordinals_b.iter().chain(predicted_fragment_ordinals_y.iter()).cloned().collect(),
            intensities: predicted_intensities_b.iter().chain(predicted_intensities_y.iter()).cloned().collect(),
            mz_calculated: fragments_observed.inner.mz_calculated.clone(),
            mz_experimental: fragments_observed.inner.mz_experimental.clone(),
        }
    };

    let fragment_intensity_pred = FragmentIntensityPrediction::new(
        fragments_observed.inner.intensities.clone(),
        fragments_observed.inner.mz_experimental.clone(),
        fragments_observed.inner.mz_calculated.clone(),
        fragments_observed.inner.charges.clone(),
        fragments_observed.inner.fragment_ordinals.clone(),
        fragments_observed.inner.kinds.iter().map(|k| kind_to_string(*k) == "y").collect(),
        intensity_copy.clone(),
    );

    let cosine_sim = fragment_intensity_pred.cosine_similarity(1e-7, false);
    let spearman = fragment_intensity_pred.spearman_correlation(1e-7, false);
    let pearson = fragment_intensity_pred.pearson_correlation(1e-7, false);
    let spectral_entropy = fragment_intensity_pred.spectral_entropy_similarity(1e-7, false);

    let mut psm = PyPeptideSpectrumMatch {
        inner: psm.inner,
        fragments_observed: Some(fragments_observed.clone()),
        fragments_predicted: Some(fragments_predicted),
    };

    psm.inner.prosit_intensities = Some(intensity_copy);

    psm.inner.cosine_similarity = cosine_sim;
    psm.inner.spectral_entropy_similarity = Some(spectral_entropy);
    psm.inner.spectral_correlation_similarity_pearson = Some(pearson);
    psm.inner.spectral_correlation_similarity_spearman = Some(spearman);

    psm
}

#[pyfunction]
pub fn associate_fragment_ions_with_prosit_predicted_intensities_par(
    psms: Vec<PyPeptideSpectrumMatch>,
    flat_intensities: Vec<Vec<f32>>,
    num_threads: usize
) -> Vec<PyPeptideSpectrumMatch> {
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let result = pool.install(|| {
        psms.par_iter().zip(flat_intensities.par_iter())
            .map(|(psm, flat_intensities)| {
                associate_psm_with_prosit_predicted_intensities(psm.clone(), flat_intensities.clone())
            })
            .collect()
    });

    result
}

#[pyfunction]
pub fn merge_psm_maps(left_map: BTreeMap<String, Vec<PyPsm>>, right_map: BTreeMap<String, Vec<PyPsm>>, max_hits: usize) -> BTreeMap<String, Vec<PyPsm>> {

    // generate correct peptide status and protein ids
    let peptide_map = get_peptide_map(left_map.clone(), right_map.clone());

    // update left and right map with decoy status and protein ids from the peptide map
    let left_map = update_psm_map(left_map, peptide_map.clone());
    let right_map = update_psm_map(right_map, peptide_map);

    // merge the two maps
    let mut merged_map: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();

    for (key, psms) in left_map.iter().chain(right_map.iter()) {
        if !merged_map.contains_key(key) {
            merged_map.insert(key.clone(), Vec::new());
        }
        merged_map.get_mut(key).unwrap().extend(psms.clone());
    }

    // remove duplicates
    merged_map = remove_duplicates(merged_map);

    // clip the number of hits
    for (_, psms) in merged_map.iter_mut() {
        psms.truncate(max_hits);
    }

    merged_map
}

fn update_psm_map(psm_map: BTreeMap<String, Vec<PyPsm>>, peptide_map: BTreeMap<String, (bool, Vec<String>)>) -> BTreeMap<String, Vec<PyPsm>> {
    let mut new_map: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();

    for (key, psms) in psm_map {
        let mut new_psms: Vec<PyPsm> = Vec::new();
        for psm in psms {
            let sequence = psm.clone().inner.sequence.unwrap().sequence;

            // if the peptide is not in the map, skip the psm
            if peptide_map.get(&sequence).is_none() {
                continue;
            }

            let (decoy, proteins) = peptide_map.get(&sequence).unwrap();
            let mut new_psm = psm.clone();
            new_psm.inner.decoy = *decoy;
            new_psm.inner.proteins = proteins.clone();
            new_psms.push(new_psm);
        }
        new_map.insert(key, new_psms);
    }

    new_map
}

fn remove_duplicates(psm_map: BTreeMap<String, Vec<PyPsm>>) -> BTreeMap<String, Vec<PyPsm>> {
    let mut new_map: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();

    for (key, psms) in psm_map {
        let mut new_psms: Vec<PyPsm> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        // sort the psms by hyperscore descending
        for psm in psms.iter().sorted_by(|a, b| b.inner.hyperscore.partial_cmp(&a.inner.hyperscore).unwrap()) {
            let sequence = psm.clone().inner.sequence.unwrap().sequence;
            if !seen.contains(&sequence) {
                seen.insert(sequence.clone());
                new_psms.push(psm.clone());
            }
        }
        new_map.insert(key, new_psms);
    }

    new_map
}

fn get_peptide_map(left_map: BTreeMap<String, Vec<PyPsm>>, right_map: BTreeMap<String, Vec<PyPsm>>) -> BTreeMap<String, (bool, Vec<String>)> {

    let mut peptide_map: BTreeMap<String, (bool, HashSet<String>)> = BTreeMap::new();

    let psms = left_map.into_iter().chain(right_map.into_iter()).collect::<Vec<_>>();

    for (_, psms) in psms {
        for psm in psms {
            let key = psm.inner.sequence.unwrap().sequence;
            let decoy = psm.inner.decoy;
            let proteins = psm.inner.proteins;

            // if the peptide is already in the map
            if peptide_map.contains_key(&key) {

                let (current_decoy, current_proteins) = peptide_map.get_mut(&key).unwrap();

                // if decoy of the current peptide is false and the new decoy is also false, add the proteins to the current proteins
                if !decoy && !*current_decoy {
                    current_proteins.extend(proteins);
                }

                // if both are decoy, merge the proteins
                else if decoy && *current_decoy {
                    current_proteins.extend(proteins);
                }

                // if new is not decoy, set the peptide to the new peptide and the proteins to the new proteins
                else if !decoy && *current_decoy {
                    *current_decoy = decoy;
                    *current_proteins = proteins.into_iter().collect();
                }

                // if the new is decoy but the current is not decoy, do nothing
                else if decoy && !*current_decoy {
                    continue;
                }

                /*
                // if new is not decoy but current is decoy, remove this peptide from the map since no overlap between decoy and target should exist
                else if !decoy && *current_decoy {
                    peptide_map.remove(&key);
                }

                // also, if the peptide is not decoy but the current peptide is decoy, remove the current peptide from the map
                else if decoy && !*current_decoy {
                    peptide_map.remove(&key);
                }
                 */

            // if the peptide is not in the map, add it
            } else {
                peptide_map.insert(key, (decoy, proteins.into_iter().collect()));
            }
        }
    }

    // return the peptide map by converting the hashset to a vector
    peptide_map.into_iter().map(|(k, (d, p))| (k, (d, p.into_iter().collect()))).collect()
}

fn kind_to_string(kind: Kind) -> String {
    match kind {
        Kind::B => "b".to_string(),
        Kind::Y => "y".to_string(),
        Kind::Z => "z".to_string(),
        Kind::A => "a".to_string(),
        Kind::C => "c".to_string(),
        Kind::X => "x".to_string(),
    }
}

fn string_to_kind(kind: &str) -> Kind {
    match kind {
        "b" => Kind::B,
        "y" => Kind::Y,
        "z" => Kind::Z,
        "a" => Kind::A,
        "c" => Kind::C,
        "x" => Kind::X,
        _ => panic!("Invalid kind"),
    }
}

#[pymodule]
pub fn scoring(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragments>()?;
    m.add_class::<PyFeature>()?;
    m.add_class::<PyScorer>()?;
    m.add_class::<PyPsm>()?;
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_class::<PyScoreType>()?;
    m.add_function(wrap_pyfunction!(associate_psm_with_prosit_predicted_intensities, m)?)?;
    m.add_function(wrap_pyfunction!(associate_fragment_ions_with_prosit_predicted_intensities_par, m)?)?;
    m.add_function(wrap_pyfunction!(prosit_intensities_to_py_fragments, m)?)?;
    m.add_function(wrap_pyfunction!(prosit_intensities_to_py_fragments_par, m)?)?;
    m.add_function(wrap_pyfunction!(psm_from_json, m)?)?;
    m.add_function(wrap_pyfunction!(merge_psm_maps, m)?)?;
    m.add_function(wrap_pyfunction!(peptide_spectrum_match_from_json, m)?)?;
    Ok(())
}