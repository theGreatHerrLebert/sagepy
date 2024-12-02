use std::collections::HashMap;
use pyo3::prelude::*;
use pyo3::types::PyList;
use qfdrust::dataset::{TDCMethod};
use qfdrust::picked::{Row};
use qfdrust::psm::Psm;
use crate::py_scoring::PyPsm;

#[pyclass]
#[derive(Clone)]
pub struct PyRow {
    pub inner: Row
}

#[pymethods]
impl PyRow {
    #[new]
    fn new(spec_idx: String, match_idx: String, decoy: bool, score: f32, q_value: f64) -> Self {
        PyRow {
            inner: Row {
                key: (spec_idx, match_idx),
                decoy,
                score,
                q_value
            }
        }
    }

    #[getter]
    fn spec_idx(&self) -> &str {
        &self.inner.key.0
    }

    #[getter]
    fn match_idx(&self) -> &str {
        &self.inner.key.1
    }

    #[getter]
    fn decoy(&self) -> bool {
        self.inner.decoy
    }

    #[getter]
    fn score(&self) -> f32 {
        self.inner.score
    }
}

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

#[pyfunction]
pub fn target_decoy_competition(
    method: &PyTDCMethod,
    spectra_idx: Vec<String>,
    match_idx: Vec<String>,
    target: Vec<bool>,
    scores: Vec<f32>,
    match_identiy_candidates: Vec<Option<Vec<String>>>,
) -> (Vec<String>, Vec<String>, Vec<Vec<String>>, Vec<bool>, Vec<f32>, Vec<f64>) {

    let method = method.inner.clone();

    let (spec_idx, match_idx, match_identity, decoy, scores, q_values) = qfdrust::dataset::target_decoy_competition(method, spectra_idx, match_idx, target, scores, match_identiy_candidates);

    (spec_idx, match_idx, match_identity, decoy, scores, q_values)
}

#[pyfunction]
pub fn assign_spectrum_q(_py: Python, psm_collection: &PyList, use_hyper_score: bool) -> PyResult<()> {
    // Extract the inner collection of Feature objects along with their original indices
    let mut inner_collection: Vec<Psm> = psm_collection.iter().map(|item| {
            // Extract each item as a PyCell<PyPsm>
            let feature: &PyCell<PyPsm> = item.extract().expect("Failed to extract PyPsm");
            // Clone the inner Feature and keep the original index
            feature.borrow().inner.clone()
        }).collect();

    let q_values = qfdrust::picked::spectrum_q_value(&inner_collection, use_hyper_score);

    // Update the q_values in the inner collection
    for (psm, q_value) in inner_collection.iter_mut().zip(q_values.iter()) {
        psm.sage_feature.spectrum_q = *q_value;
    }

    Ok(())
}

#[pyfunction]
pub fn assign_q_values(
    rows: Vec<PyRow>,
) -> HashMap<(String, String), f64> {
    let rows = rows.iter().map(|row| row.inner.clone()).collect();
    qfdrust::picked::assign_q_value(rows)
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTDCMethod>()?;
    m.add_function(wrap_pyfunction!(target_decoy_competition, m)?)?;
    m.add_function(wrap_pyfunction!(assign_q_values, m)?)?;
    m.add_function(wrap_pyfunction!(assign_spectrum_q, m)?)?;
    m.add_class::<PyRow>()?;
    Ok(())
}