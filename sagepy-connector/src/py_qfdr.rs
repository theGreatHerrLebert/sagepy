use pyo3::prelude::*;
use qfdrust::dataset::{TDCMethod};

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
    scores: Vec<f32>) -> (Vec<String>, Vec<String>, Vec<bool>, Vec<f32>, Vec<f64>) {

    let method = method.inner.clone();

    let (spec_idx, match_idx, decoy, scores, q_values) = qfdrust::dataset::target_decoy_competition(method, spectra_idx, match_idx, target, scores);

    (spec_idx, match_idx, decoy, scores, q_values)
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTDCMethod>()?;
    m.add_function(wrap_pyfunction!(target_decoy_competition, m)?)?;
    Ok(())
}