use pyo3::prelude::*;
use qfdrust::dataset::{PeptideSpectrumMatch, TDCMethod};
use sage_core::scoring::Fragments;
use crate::py_scoring::PyFragments;

#[pyclass]
#[derive(Clone)]
pub struct PyTDCMethod {
    pub inner: TDCMethod,
}

#[pymethods]
impl PyTDCMethod {
    #[new]
    fn new(method: &str) -> Self {
        PyTDCMethod {
            inner: TDCMethod::from_str(method),
        }
    }
    pub fn to_str(&self) -> &str {
        self.inner.to_str()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyPeptideSpectrumMatch {
    pub inner: PeptideSpectrumMatch,
    pub fragments: Option<Fragments>,
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
        sequence: Option<String>,
        charge: Option<u8>,
        retention_time_observed: Option<f32>,
        retention_time_predicted: Option<f32>,
        inverse_mobility_observed: Option<f32>,
        inverse_mobility_predicted: Option<f32>,
        intensity_ms1: Option<f32>,
        intensity_ms2: Option<f32>,
        q_value: Option<f64>,
        fragments: Option<PyFragments>,
    ) -> Self {

        let fragments = match fragments {
            Some(fragments) => Some(fragments.inner),
            None => None,
        };

        PyPeptideSpectrumMatch {
            inner: PeptideSpectrumMatch::new(
                spec_idx,
                peptide_idx,
                proteins,
                decoy,
                hyper_score,
                rank,
                mono_mass_observed,
                sequence,
                charge,
                retention_time_observed,
                retention_time_predicted,
                inverse_mobility_observed,
                inverse_mobility_predicted,
                intensity_ms1,
                intensity_ms2,
                q_value,
            ),
            fragments,
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

    #[getter]
    pub fn rank(&self) -> u32 {
        self.inner.rank
    }

    #[setter]
    pub fn set_hyper_score(&mut self, hyper_score: f64) {
        self.inner.hyper_score = hyper_score;
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
    pub fn fragments_predicted(&self) -> Option<PyFragments> {
        todo!("Implement this")
    }

    #[getter]
    pub fn fragments_observed(&self) -> Option<PyFragments> {
        let maybe_fragments = self.fragments.clone();
        match maybe_fragments {
            Some(fragments) => Some(PyFragments { inner: fragments }),
            None => None,
        }
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

    pub fn associate_fragment_ions_with_prosit_predicted_intensities(&mut self, flat_intensities: Vec<f64>) {
        let ion_series = self.inner.associate_with_prosit_predicted_intensities(flat_intensities);
        self.inner.peptide_product_ion_series_collection_predicted = ion_series;
    }
}

#[pyfunction]
pub fn target_decoy_competition(
    method: &PyTDCMethod,
    spectra_idx: Vec<String>,
    match_idx: Vec<u32>,
    target: Vec<bool>,
    scores: Vec<f32>) -> (Vec<String>, Vec<u32>, Vec<bool>, Vec<f32>, Vec<f64>) {

    let method = method.inner.clone();

    let (spec_idx, match_idx, decoy, scores, q_values) = qfdrust::dataset::target_decoy_competition(method, spectra_idx, match_idx, target, scores);

    (spec_idx, match_idx, decoy, scores, q_values)
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_class::<PyTDCMethod>()?;
    m.add_function(wrap_pyfunction!(target_decoy_competition, m)?)?;
    Ok(())
}