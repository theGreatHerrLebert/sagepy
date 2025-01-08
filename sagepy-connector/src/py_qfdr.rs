use pyo3::prelude::*;
use pyo3::types::PyList;
use qfdrust::dataset::{TDCMethod};
use qfdrust::picked::protein_id_from_psm;
use qfdrust::psm::Psm;
use crate::py_scoring::PyPsm;

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
    let inner_collection: Vec<Psm> = psm_collection.iter().map(|item| {
            // Extract each item as a PyCell<PyPsm>
            let feature: &PyCell<PyPsm> = item.extract().expect("Failed to extract PyPsm");
            // Clone the inner Feature and keep the original index
            feature.borrow().inner.clone()
        }).collect();

    let q_values = qfdrust::picked::spectrum_q_value(&inner_collection, use_hyper_score);

    // Update the q_values of the inner collection IN PLACE, meaning its necessarry to fetch them by index for PyO3 to work
    for (index, q_value) in q_values.iter().enumerate() {
        let feature: &PyCell<PyPsm> = psm_collection.get_item(index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();
        feature_borrow.inner.sage_feature.spectrum_q = q_value.clone();
    }

    Ok(())
}

#[pyfunction]
pub fn assign_peptide_q(_py: Python, psm_collection: &PyList, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut  inner_collection: Vec<Psm> = psm_collection.iter().map(|item| {
            // Extract each item as a PyCell<PyPsm>
            let feature: &PyCell<PyPsm> = item.extract().expect("Failed to extract PyPsm");
            // Clone the inner Feature and keep the original index
            feature.borrow().inner.clone()
        }).collect();

    let q_values = qfdrust::picked::picked_peptide(&mut inner_collection, use_hyper_score);

    // Update the q_values of the inner collection IN PLACE, meaning its necessary to fetch them by index for PyO3 to work
    for (index, _) in psm_collection.iter().enumerate() {
        let feature: &PyCell<PyPsm> = psm_collection.get_item(index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();

        let key = match feature_borrow.inner.sage_feature.label {
            -1 => feature_borrow.inner.sequence_decoy.clone().unwrap().sequence.clone(),
            _ => feature_borrow.inner.sequence.clone().unwrap().sequence.clone(),
        };

        feature_borrow.inner.sage_feature.peptide_q = q_values.get(&key).unwrap_or(&1.0).clone() as f32;
    }

    Ok(())
}

#[pyfunction]
pub fn assign_protein_q(_py: Python, psm_collection: &PyList, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut inner_collection: Vec<Psm> = psm_collection.iter().map(|item| {
            // Extract each item as a PyCell<PyPsm>
            let feature: &PyCell<PyPsm> = item.extract().expect("Failed to extract PyPsm");
            // Clone the inner Feature and keep the original index
            feature.borrow().inner.clone()
        }).collect();

    let q_values = qfdrust::picked::picked_protein(&mut inner_collection, use_hyper_score);

    // Update the q_values of the inner collection IN PLACE, meaning its necessarry to fetch them by index for PyO3 to work
    for (index, _) in psm_collection.iter().enumerate() {
        let feature: &PyCell<PyPsm> = psm_collection.get_item(index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();

        let key = protein_id_from_psm(&feature_borrow.inner, "rev_", true);

        feature_borrow.inner.sage_feature.protein_q = q_values.get(&key).unwrap_or(&1.0).clone() as f32;
    }

    Ok(())
}

#[pymodule]
pub fn py_qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTDCMethod>()?;
    m.add_function(wrap_pyfunction!(target_decoy_competition, m)?)?;
    m.add_function(wrap_pyfunction!(assign_spectrum_q, m)?)?;
    m.add_function(wrap_pyfunction!(assign_peptide_q, m)?)?;
    m.add_function(wrap_pyfunction!(assign_protein_q, m)?)?;
    Ok(())
}