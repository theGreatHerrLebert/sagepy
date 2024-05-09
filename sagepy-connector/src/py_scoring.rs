use std::collections::{BTreeMap, HashSet};
use pyo3::prelude::*;
use qfdrust::dataset::{PeptideSpectrumMatch};
use qfdrust::utility::sage_sequence_to_unimod_sequence;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use sage_core::ion_series::Kind;

use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_mass::PyTolerance;
use crate::py_spectrum::{PyProcessedSpectrum};
use sage_core::scoring::{Feature, Scorer, Fragments};
use crate::py_ion_series::PyKind;
use crate::py_utility::{cosine_similarity, flat_prosit_array_to_fragments_map, py_fragments_to_fragments_map};

#[pyclass]
#[derive(Clone)]
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
    pub min_fragment_mass: f32,
    pub max_fragment_mass: f32,
    pub chimera: bool,
    pub report_psms: usize,
    pub wide_window: bool,
    pub annotate_matches: bool,
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
        min_fragment_mass: f32,
        max_fragment_mass: f32,
        chimera: bool,
        report_psms: usize,
        wide_window: bool,
        annotate_matches: bool,
        max_fragment_charge: Option<u8>,
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
            min_fragment_mass,
            max_fragment_mass,
            chimera,
            report_psms,
            wide_window,
            annotate_matches,
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
            min_fragment_mass: self.min_fragment_mass,
            max_fragment_mass: self.max_fragment_mass,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
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
            min_fragment_mass: self.min_fragment_mass,
            max_fragment_mass: self.max_fragment_mass,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
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

    pub fn score_collection_to_psm_collection(
        &self,
        db: &PyIndexedDatabase,
        spectra: Vec<PyProcessedSpectrum>,
        num_threads: usize,
    ) -> Vec<PyPeptideSpectrumMatch> {
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
            min_fragment_mass: self.min_fragment_mass,
            max_fragment_mass: self.max_fragment_mass,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
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
                        let intensity_ms1: f32 = spectrum.inner.precursors.iter().map(|p| p.intensity.unwrap()).sum();
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
                            Some(sage_sequence_to_unimod_sequence(sequence, &peptide.modifications)),
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

        // sort by sec_id, match_id, decoy
        result.sort_by(|a, b| {
            let a = &a.inner;
            let b = &b.inner;
            a.spec_idx.cmp(&b.spec_idx)
                .then(a.peptide_idx.cmp(&b.peptide_idx))
                .then(a.decoy.cmp(&b.decoy))
        });

        result
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
            min_fragment_mass: self.min_fragment_mass,
            max_fragment_mass: self.max_fragment_mass,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
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
            min_fragment_mass: self.min_fragment_mass,
            max_fragment_mass: self.max_fragment_mass,
            chimera: self.chimera,
            report_psms: self.report_psms,
            wide_window: self.wide_window,
            annotate_matches: self.annotate_matches,
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
    pub fn min_fragment_mass(&self) -> f32 {
        self.min_fragment_mass
    }

    #[getter]
    pub fn max_fragment_mass(&self) -> f32 {
        self.max_fragment_mass
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
#[derive(Clone)]
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
    ) -> Self {

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
    pub fn fragments_observed(&self) -> Option<PyFragments> {
        self.fragments_observed.clone()
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
        match (&self.fragments_observed, &self.fragments_predicted) {
            (Some(observed), Some(predicted)) => {
                let observed_intensities = observed.intensities();
                let predicted_intensities = predicted.intensities();
                cosine_similarity(&observed_intensities, &predicted_intensities)
            }
            _ => None,
        }
    }
}

#[pyfunction]
pub fn associate_psm_with_prosit_predicted_intensities(
    psm: PyPeptideSpectrumMatch,
    flat_intensities: Vec<f32>,
) -> PyPeptideSpectrumMatch {

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

    PyPeptideSpectrumMatch {
        inner: psm.inner,
        fragments_observed: Some(fragments_observed.clone()),
        fragments_predicted: Some(fragments_predicted),
    }
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
        psms.par_iter()
            .zip(flat_intensities.par_iter())
            .map(|(psm, flat_intensities)| {
                associate_psm_with_prosit_predicted_intensities(psm.clone(), flat_intensities.clone())
            })
            .collect()
    });

    result
}

#[pymodule]
pub fn scoring(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragments>()?;
    m.add_class::<PyFeature>()?;
    m.add_class::<PyScorer>()?;
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_function(wrap_pyfunction!(associate_psm_with_prosit_predicted_intensities, m)?)?;
    m.add_function(wrap_pyfunction!(associate_fragment_ions_with_prosit_predicted_intensities_par, m)?)?;
    Ok(())
}