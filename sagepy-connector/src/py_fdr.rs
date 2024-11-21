use pyo3::prelude::*;
use sage_core::fdr::{Competition, picked_peptide, picked_protein};
use sage_core::database::{PeptideIx};
use sage_core::scoring::Feature;
use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_scoring::{PyFeature, PyPsm};

#[pyclass]
// TODO: Check if it makes sense to tie this to PeptideIx
struct PyCompetitionPeptideIx {
    inner: Competition<PeptideIx>,
}

#[pymethods]
impl PyCompetitionPeptideIx {
    #[new]
    fn new(forward: f32, reverse: f32, forward_ix: Option<PyPeptideIx>, reverse_ix: Option<PyPeptideIx>) -> Self {
        PyCompetitionPeptideIx {
            inner: Competition {
                forward,
                foward_ix: forward_ix.map(|ix| ix.inner),
                reverse,
                reverse_ix: reverse_ix.map(|ix| ix.inner),
            },
        }
    }
    #[getter]
    fn forward(&self) -> f32 {
        self.inner.forward
    }

    #[getter]
    fn reverse(&self) -> f32 {
        self.inner.reverse
    }

    #[getter]
    fn forward_ix(&self) -> Option<PyPeptideIx> {
        self.inner.foward_ix.map(|ix| PyPeptideIx { inner: ix })
    }

    #[getter]
    fn reverse_ix(&self) -> Option<PyPeptideIx> {
        self.inner.reverse_ix.map(|ix| PyPeptideIx { inner: ix })
    }
}

#[pyfunction]
pub fn py_picked_peptide(mut feature_collection: Vec<PyFeature>, indexed_database: &PyIndexedDatabase) {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.clone()).collect();
    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.peptide_q = inner.peptide_q;
        feature.inner.protein_q = inner.protein_q;
        feature.inner.posterior_error = inner.posterior_error;
    }
}

#[pyfunction]
pub fn py_picked_protein(mut feature_collection: Vec<PyFeature>, indexed_database: &PyIndexedDatabase) {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.clone()).collect();
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.peptide_q = inner.peptide_q;
        feature.inner.protein_q = inner.protein_q;
        feature.inner.posterior_error = inner.posterior_error;
    }
}

#[pyfunction]
pub fn py_picked_peptide_psm(mut feature_collection: Vec<PyPsm>, indexed_database: &PyIndexedDatabase) {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.sage_feature.clone()).collect();
    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.sage_feature.peptide_q = inner.peptide_q;
        feature.inner.sage_feature.protein_q = inner.protein_q;
        feature.inner.sage_feature.posterior_error = inner.posterior_error;
    }
}

#[pyfunction]
pub fn py_picked_protein_psm(mut feature_collection: Vec<PyPsm>, indexed_database: &PyIndexedDatabase) {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.sage_feature.clone()).collect();
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.sage_feature.peptide_q = inner.peptide_q;
        feature.inner.sage_feature.protein_q = inner.protein_q;
        feature.inner.sage_feature.posterior_error = inner.posterior_error;
    }
}

#[pymodule]
pub fn fdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    m.add_function(wrap_pyfunction!(py_picked_peptide, m)?)?;
    m.add_function(wrap_pyfunction!(py_picked_protein, m)?)?;
    m.add_function(wrap_pyfunction!(py_picked_peptide_psm, m)?)?;
    m.add_function(wrap_pyfunction!(py_picked_protein_psm, m)?)?;
    Ok(())
}