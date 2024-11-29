use std::collections::HashMap;
use pyo3::prelude::*;
use qfdrust::dataset::{TDCMethod};
use qfdrust::picked::{Row};

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
    m.add_class::<PyRow>()?;
    Ok(())
}