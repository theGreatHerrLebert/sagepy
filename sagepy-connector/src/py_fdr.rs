use pyo3::prelude::*;
use sage_core::fdr::{Competition, picked_peptide, picked_protein};
use sage_core::database::{IndexedDatabase, PeptideIx};
use sage_core::scoring::Feature;
use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_scoring::PyFeature;

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
pub fn py_picked_peptide(feature_collection: Vec<PyFeature>, indexed_database: &PyIndexedDatabase) -> Vec<PyFeature> {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.clone()).collect();
    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);
    inner_collection.into_iter().map(|feature| PyFeature { inner: feature }).collect()
}

#[pyfunction]
pub fn py_picked_protein(feature_collection: Vec<PyFeature>, indexed_database: &PyIndexedDatabase) -> Vec<PyFeature> {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.clone()).collect();
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);
    inner_collection.into_iter().map(|feature| PyFeature { inner: feature }).collect()
}

#[pymodule]
pub fn fdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    m.add_function(wrap_pyfunction!(py_picked_peptide, m)?)?;
    m.add_function(wrap_pyfunction!(py_picked_protein, m)?)?;
    Ok(())
}