use numpy::inner;
use pyo3::prelude::*;
use sage_core::fdr::{Competition, picked_peptide, picked_protein};
use sage_core::database::{PeptideIx};
use sage_core::scoring::Feature;
use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_scoring::PyFeature;
use rayon::prelude::*;

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

    inner_collection.par_iter_mut().for_each(|feat| {
        feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
    });

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

    inner_collection.par_iter_mut().for_each(|feat| {
        feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
    });

    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.peptide_q = inner.peptide_q;
        feature.inner.protein_q = inner.protein_q;
        feature.inner.posterior_error = inner.posterior_error;
    }
}
#[pyfunction]
pub fn py_sage_fdr(mut feature_collection: Vec<PyFeature>, indexed_database: &PyIndexedDatabase) {
    let mut inner_collection: Vec<Feature> = feature_collection.iter().map(|feature| feature.inner.clone()).collect();

    inner_collection.par_iter_mut().for_each(|feat| {
        feat.discriminant_score = feat.hyperscore as f32//(-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
    });

    inner_collection.par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));

    let inner_collection_mut = inner_collection.as_mut_slice();

    let _ = sage_core::ml::qvalue::spectrum_q_value(inner_collection_mut);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection_mut.iter()) {
        feature.inner.discriminant_score = inner.discriminant_score;
        feature.inner.spectrum_q = inner.spectrum_q;
    }

    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    for (feature, inner) in feature_collection.iter_mut().zip(inner_collection.iter()) {
        feature.inner.spectrum_q = inner.spectrum_q;
        feature.inner.peptide_q = inner.peptide_q;
        feature.inner.protein_q = inner.protein_q;
        feature.inner.posterior_error = inner.posterior_error;
    }
}

#[pymodule]
pub fn fdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    m.add_function(wrap_pyfunction!(py_picked_peptide, m)?)?;
    m.add_function(wrap_pyfunction!(py_picked_protein, m)?)?;
    m.add_function(wrap_pyfunction!(py_sage_fdr, m)?)?;
    Ok(())
}