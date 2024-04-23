use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::dataset::{PeptideSpectrumMatch, PsmDataset, TDCMethod};

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
        mono_mass_calculated: Option<f32>,
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
        re_score: Option<f64>,
    ) -> Self {
        PyPeptideSpectrumMatch {
            inner: PeptideSpectrumMatch::new(
                spec_idx,
                peptide_idx,
                proteins,
                decoy,
                hyper_score,
                mono_mass_calculated,
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
                re_score,
            ),
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
}

#[pyclass]
pub struct PyPsmDataset {
    pub inner: PsmDataset,
}

#[pymethods]
impl PyPsmDataset {
    #[new]
    fn new(spec_ids: Vec<String>, matches: Vec<Vec<PyPeptideSpectrumMatch>>) -> Self {
        let mut psm_map = BTreeMap::new();
        let inner_matches: Vec<Vec<_>> = matches.into_iter().map(|m| m.iter().map(|m| m.inner.clone()).collect()).collect();
        for (spec_id, matches) in spec_ids.into_iter().zip(inner_matches.into_iter()) {
            psm_map.insert(spec_id, matches);
        }
        PyPsmDataset {
            inner: PsmDataset {
                psm_map,
            },
        }
    }
    pub fn get_spec_psms(&self, spec_id: String) -> PyResult<Vec<PyPeptideSpectrumMatch>> {
        Ok(self.inner.psm_map.get(&spec_id).unwrap().iter().map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() }).collect())
    }

    #[getter]
    pub fn size(&self) -> usize {
        self.inner.size()
    }

    #[getter]
    pub fn keys(&self) -> Vec<String> {
        self.inner.get_spectra_ids()
    }

    pub fn best_target_ptm(&self, spec_id: String) -> Option<PyPeptideSpectrumMatch> {
        self.inner.get_best_target_psm(&spec_id).map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() })
    }

    pub fn best_decoy_ptm(&self, spec_id: String) -> Option<PyPeptideSpectrumMatch> {
        self.inner.get_best_decoy_psm(&spec_id).map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() })
    }

    pub fn tdc(&self, method: PyTDCMethod) -> Vec<PyPeptideSpectrumMatch> {
        self.inner.tdc(method.inner).iter().map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() }).collect()
    }
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_class::<PyPsmDataset>()?;
    m.add_class::<PyTDCMethod>()?;
    Ok(())
}