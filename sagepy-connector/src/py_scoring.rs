use std::collections::{BTreeMap, HashSet};
use itertools::Itertools;
use pyo3::prelude::*;
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
    #[pyo3(signature = (spec_idx, peptide_idx, proteins, sage_feature, sequence=None, sequence_modified=None, sequence_decoy=None, sequence_decoy_modified=None, intensity_ms1=None, intensity_ms2=None, collision_energy=None, collision_energy_calibrated=None, retention_time_projected=None, prosit_predicted_intensities=None, re_score=None))]
    pub fn new(
        spec_idx: String,
        peptide_idx: u32,
        proteins: Vec<String>,
        sage_feature: PyFeature,
        sequence: Option<String>,
        sequence_modified: Option<String>,
        sequence_decoy: Option<String>,
        sequence_decoy_modified: Option<String>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        collision_energy: Option<f32>,
        collision_energy_calibrated: Option<f32>,
        retention_time_projected: Option<f32>,
        prosit_predicted_intensities: Option<Vec<f32>>,
        re_score: Option<f64>,
    ) -> Self {
        PyPsm {
            inner: Psm::new(
                spec_idx,
                peptide_idx,
                proteins,
                sage_feature.inner.clone(),
                sequence,
                sequence_modified,
                sequence_decoy,
                sequence_decoy_modified,
                intensity_ms1,
                intensity_ms2,
                collision_energy,
                collision_energy_calibrated,
                retention_time_projected,
                prosit_predicted_intensities,
                re_score,
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
        self.inner.sage_feature.hyperscore
    }
    
    #[setter]
    pub fn set_hyperscore(&mut self, value: f64) {
        self.inner.sage_feature.hyperscore = value;
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

    #[setter]
    pub fn set_sage_feature_file_id(&mut self, value: usize) {
        self.inner.sage_feature.file_id = value;
    }

    #[getter]
    pub fn sequence(&self) -> Option<String> {
        Some(self.inner.clone().sequence?.sequence)
    }

    #[getter]
    pub fn sequence_modified(&self) -> Option<String> {
        Some(self.inner.clone().sequence_modified?.sequence)
    }

    #[getter]
    pub fn sequence_decoy(&self) -> Option<String> {
        Some(self.inner.clone().sequence_decoy?.sequence)
    }

    #[getter]
    pub fn sequence_decoy_modified(&self) -> Option<String> {
        Some(self.inner.clone().sequence_decoy_modified?.sequence)
    }

    #[getter]
    pub fn charge(&self) -> u8 {
        self.inner.sage_feature.charge
    }
    
    #[setter]
    pub fn set_charge(&mut self, value:u8) {
        self.inner.sage_feature.charge = value;
    }
    
    #[getter]
    pub fn mono_mass_observed(&self) -> f32 {
        self.inner.sage_feature.expmass
    }
    
    #[setter]
    pub fn set_mono_mass_observed(&mut self, value: f32) {
        self.inner.sage_feature.expmass = value;
    }

    #[getter]
    pub fn mono_mass_calculated(&self) -> f32 {
        self.inner.sage_feature.calcmass
    }
    
    #[setter]
    pub fn set_mono_mass_calculated(&mut self, value:f32) {
        self.inner.sage_feature.calcmass = value;
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
    pub fn retention_time(&self) -> f32 {
        self.inner.sage_feature.rt
    }
    
    #[setter]
    pub fn set_retention_time(&mut self, value: f32) {
        self.inner.sage_feature.rt = value;
    }
    
    #[getter]
    pub fn retention_time_calibrated(&self) -> f32 {
        self.inner.sage_feature.aligned_rt
    }
    
    #[setter]
    pub fn set_retention_time_calibrated(&mut self, value: f32) {
        self.inner.sage_feature.aligned_rt = value;
    }

    #[getter]
    pub fn retention_time_predicted(&self) -> f32 {
        self.inner.sage_feature.predicted_rt
    }

    #[setter]
    pub fn set_retention_time_predicted(&mut self, value: f32) {
        self.inner.sage_feature.predicted_rt = value;
        self.inner.sage_feature.delta_rt_model = value - self.inner.retention_time_projected.unwrap_or(0.0);
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
    pub fn inverse_ion_mobility(&self) -> f32 {
        self.inner.sage_feature.ims
    }
    
    #[setter]
    pub fn set_inverse_ion_mobility(&mut self, value:f32) {
        self.inner.sage_feature.ims = value;
    }
    
    #[getter]
    pub fn inverse_ion_mobility_predicted(&self) -> f32 {
        self.inner.sage_feature.predicted_ims
    }
    
    #[setter]
    pub fn set_inverse_ion_mobility_predicted(&mut self, value: f32) {
        self.inner.sage_feature.predicted_ims = value;
        self.inner.sage_feature.delta_ims_model = value - self.inner.sage_feature.ims;
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
    pub fn rank(&self) -> u32 {
        self.inner.sage_feature.rank
    }

    #[setter]
    pub fn set_rank(&mut self, value: u32) {
        self.inner.sage_feature.rank = value;
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
    pub fn decoy(&self) -> bool {
        self.inner.sage_feature.label == -1
    }
    
    #[setter]
    pub fn set_decoy(&mut self, value: bool) {
        self.inner.sage_feature.label = if value { -1 } else { 1 };
    }

    #[getter]
    pub fn get_spectral_angle_similarity(&self) -> f32 {
        self.inner.fragment_intensity_prediction.clone().unwrap().spectral_angle_similarity(0.001, false)
    }

    pub fn get_fragment_intensity_prediction(&self) -> PyFragmentIntensityPrediction {
        PyFragmentIntensityPrediction {
            inner: self.inner.fragment_intensity_prediction.clone().unwrap(),
        }
    }

    pub fn get_feature_names(&self) -> Vec<&str> {
        self.inner.get_feature_names()
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(&self.inner).unwrap()
    }

    pub fn prosit_intensities_to_fragments(&self) -> PyFragments {
        let fragments = self.inner.prosit_intensity_to_fragments();
        PyFragments {
            inner: fragments.unwrap(),
        }
    }

    pub fn observed_fragments_to_fragments_map(&self, normalize: bool) -> BTreeMap<(u32, i32, i32), f32> {
        py_fragments_to_fragments_map(&self.sage_feature().fragments().unwrap(), normalize)
    }

    pub fn prosit_intensities_to_fragments_map(&self, normalize: bool) -> BTreeMap<(u32, i32, i32), f32> {
        py_fragments_to_fragments_map(&self.prosit_intensities_to_fragments(), normalize)
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
    #[pyo3(signature = (peptide_idx, psm_id, peptide_len, spec_id, file_id, rank, label, expmass, calcmass, charge, delta_mass, isotope_error, average_ppm, hyperscore, delta_next, delta_best, matched_peaks, longest_b, longest_y, longest_y_pct, missed_cleavages, matched_intensity_pct, scored_candidates, poisson, discriminant_score, posterior_error, spectrum_q, peptide_q, protein_q, ms2_intensity, rt=None, aligned_rt=None, predicted_rt=None, delta_rt_model=None, ims=None, predicted_ims=None, delta_ims_model=None, fragments=None))]
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

    #[setter]
    pub fn set_ims(&mut self, value: f32) {
        self.inner.ims = value;
    }

    #[setter]
    pub fn set_predicted_ims(&mut self, value: f32) {
        self.inner.predicted_ims = value;
    }

    #[setter]
    pub fn set_delta_ims_model(&mut self, value: f32) {
        self.inner.delta_ims_model = value;
    }

    #[setter]
    pub fn set_rt(&mut self, value: f32) {
        self.inner.rt = value;
    }

    #[setter]
    pub fn set_aligned_rt(&mut self, value: f32) {
        self.inner.aligned_rt = value;
    }

    #[setter]
    pub fn set_predicted_rt(&mut self, value: f32) {
        self.inner.predicted_rt = value;
    }

    #[setter]
    pub fn set_delta_rt_model(&mut self, value: f32) {
        self.inner.delta_rt_model = value;
    }

    #[setter]
    pub fn set_file_id(&mut self, value: usize) {
        self.inner.file_id = value;
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
    #[pyo3(signature = (precursor_tolerance, fragment_tolerance, min_matched_peaks, min_isotope_err, max_isotope_err, min_precursor_charge, max_precursor_charge, chimera, report_psms, wide_window, annotate_matches, override_precursor_charge, expected_mods, max_fragment_charge=None, score_type=None))]
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
                        
                        let peptide = db.inner[feature.peptide_idx].clone();
                        
                        let intensity_ms1: f32 = spectrum.inner.precursors.iter().map(|p| p.intensity.unwrap_or(0.0)).sum();
                        let intensity_ms2: f32 = feature.ms2_intensity;

                        let proteins: Vec<String> = peptide.proteins.iter().map(|arc| arc.as_ref().to_string()).collect();

                        let sequence = std::str::from_utf8(&peptide.sequence).unwrap().to_string();
                        let sequence_with_mods = sage_sequence_to_unimod_sequence(sequence.clone(), &peptide.modifications, &self.expected_mods);

                        let sequence_decoy = std::str::from_utf8(&peptide.reverse(true).sequence).unwrap().to_string();
                        let sequence_decoy_with_mods = sage_sequence_to_unimod_sequence(sequence_decoy.clone(), &peptide.reverse(true).modifications, &self.expected_mods);
                        
                        let collision_energy = spectrum.collision_energies.first().unwrap_or(&None).unwrap_or(0.0f32);
                        
                        let psm = Psm::new(
                            spectrum.inner.id.clone(),
                            feature.clone().peptide_idx.0,
                            proteins,
                            feature.clone(),
                            Some(sequence), // sequence
                            Some(sequence_with_mods), // sequence_modified
                            Some(sequence_decoy), // sequence_decoy
                            Some(sequence_decoy_with_mods), // sequence_decoy_modified
                            Some(intensity_ms1),
                            Some(intensity_ms2),
                            Some(collision_energy),
                            None, // collision_energy_calibrated
                            None, // rt projected
                            None, // prosit_predicted_intensities
                            None, // re_score
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
                a.sage_feature.hyperscore.partial_cmp(&b.sage_feature.hyperscore).unwrap()
                    .then_with(|| a.peptide_idx.partial_cmp(&b.peptide_idx).unwrap())
                    .then_with(|| a.sage_feature.label.partial_cmp(&b.sage_feature.label).unwrap())
            });
        }
        
        remove_duplicates(result)
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


#[pyfunction]
pub fn psm_from_json(json: &str) -> PyPsm {
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
    psm: PyPsm,
    flat_intensities: Vec<f32>,
) -> PyPsm {
    let mut psm_copy = psm.clone();
    psm_copy.set_prosit_predicted_intensities(Some(flat_intensities.clone()));

    psm_copy
}

#[pyfunction]
pub fn associate_fragment_ions_with_prosit_predicted_intensities_par(
    psms: Vec<PyPsm>,
    flat_intensities: Vec<Vec<f32>>,
    num_threads: usize
) -> Vec<PyPsm> {
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

    let mut merged_map: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();

    // merge the maps
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

fn remove_duplicates(psm_map: BTreeMap<String, Vec<PyPsm>>) -> BTreeMap<String, Vec<PyPsm>> {
    // Parallelize the processing of the BTreeMap entries
    let new_map: BTreeMap<String, Vec<PyPsm>> = psm_map
        .into_par_iter() // Parallel iterator over the map entries
        .map(|(key, psms)| {
            let mut new_psms: Vec<PyPsm> = Vec::new();
            let mut target_seen: HashSet<String> = HashSet::new();
            let mut decoy_seen: HashSet<String> = HashSet::new();

            // Sort the psms by hyperscore descending
            for psm in psms.iter().sorted_by(|a, b| {
                b.inner
                    .sage_feature
                    .hyperscore
                    .partial_cmp(&a.inner.sage_feature.hyperscore)
                    .unwrap()
            }) {
                // Get either the target or decoy sequence
                let sequence = match psm.inner.sage_feature.label == -1 {
                    true => psm.inner.sequence_decoy.clone().unwrap().sequence,
                    false => psm.inner.sequence.clone().unwrap().sequence,
                };

                if psm.inner.sage_feature.label == -1 {
                    // Process decoys
                    if decoy_seen.contains(&sequence) {
                        let existing_psm = new_psms
                            .iter()
                            .find(|p| {
                                p.inner.sequence_decoy.clone().unwrap().sequence == sequence
                            })
                            .unwrap();
                        let mut existing_proteins = existing_psm.inner.proteins.clone();
                        for protein in &psm.inner.proteins {
                            if !existing_proteins.contains(protein) {
                                existing_proteins.push(protein.clone());
                            }
                        }
                    } else {
                        decoy_seen.insert(sequence.clone());
                        new_psms.push(psm.clone());
                    }
                } else {
                    // Process targets
                    if target_seen.contains(&sequence) {
                        let existing_psm = new_psms
                            .iter()
                            .find(|p| p.inner.sequence.clone().unwrap().sequence == sequence)
                            .unwrap();
                        let mut existing_proteins = existing_psm.inner.proteins.clone();
                        for protein in &psm.inner.proteins {
                            if !existing_proteins.contains(protein) {
                                existing_proteins.push(protein.clone());
                            }
                        }
                    } else {
                        target_seen.insert(sequence.clone());
                        new_psms.push(psm.clone());
                    }
                }
            }

            // Sort the new_psms by hyperscore descending and return the result
            let sorted_psms = new_psms
                .iter()
                .sorted_by(|a, b| {
                    b.inner
                        .sage_feature
                        .hyperscore
                        .partial_cmp(&a.inner.sage_feature.hyperscore)
                        .unwrap()
                })
                .cloned()
                .collect();

            (key, sorted_psms)
        })
        .collect();

    new_map
}
#[pyfunction]
pub fn peptide_spectrum_match_to_feature_vector(
    psm: &PyPsm,
    epsilon: f32,
    reduce_matched: bool,
) -> Vec<f32> {
    let fragment_intensity_prediction = psm.inner.get_fragment_intensity_prediction();
    fragment_intensity_prediction.get_feature_vector(epsilon, reduce_matched)
}

#[pyfunction]
pub fn peptide_spectrum_match_list_to_intensity_feature_matrix_parallel(
    psms: Vec<PyPsm>,
    epsilon: f32,
    reduce_matched: bool,
    num_threads: usize,
) -> Vec<Vec<f32>> {
    let thread_pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    thread_pool.install(|| {
        psms.par_iter().map(|psm| {
            peptide_spectrum_match_to_feature_vector(psm, epsilon, reduce_matched)
        }).collect()
    })
}

#[pymodule]
pub fn py_scoring(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFragments>()?;
    m.add_class::<PyFeature>()?;
    m.add_class::<PyScorer>()?;
    m.add_class::<PyPsm>()?;
    m.add_class::<PyScoreType>()?;
    m.add_function(wrap_pyfunction!(associate_psm_with_prosit_predicted_intensities, m)?)?;
    m.add_function(wrap_pyfunction!(associate_fragment_ions_with_prosit_predicted_intensities_par, m)?)?;
    m.add_function(wrap_pyfunction!(prosit_intensities_to_py_fragments, m)?)?;
    m.add_function(wrap_pyfunction!(prosit_intensities_to_py_fragments_par, m)?)?;
    m.add_function(wrap_pyfunction!(psm_from_json, m)?)?;
    m.add_function(wrap_pyfunction!(merge_psm_maps, m)?)?;
    m.add_function(wrap_pyfunction!(peptide_spectrum_match_to_feature_vector, m)?)?;
    m.add_function(wrap_pyfunction!(peptide_spectrum_match_list_to_intensity_feature_matrix_parallel, m)?)?;
    Ok(())
}