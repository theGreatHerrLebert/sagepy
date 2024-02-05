use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_mass::PyTolerance;
use crate::py_spectrum::PyProcessedSpectrum;
use sage_core::scoring::{Feature, Scorer, Fragments};
use crate::py_ion_series::PyKind;

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
        fragments: Option<PyFragments>
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

#[pymodule]
pub fn scoring(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragments>()?;
    m.add_class::<PyFeature>()?;
    m.add_class::<PyScorer>()?;
    Ok(())
}
