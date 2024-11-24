use numpy::inner;
use pyo3::prelude::*;
use pyo3::types::PyList;
use qfdrust::psm::Psm;
use sage_core::fdr::{Competition, picked_peptide, picked_protein};

use sage_core::database::{PeptideIx};
use sage_core::scoring::Feature;
use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_scoring::{PyFeature, PyPsm};
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
pub fn py_sage_fdr(_py: Python, feature_collection: &PyList, indexed_database: &PyIndexedDatabase, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut indexed_inner_collection: Vec<(usize, Feature)> = feature_collection.iter()
        .enumerate()
        .map(|(index, item)| {
            // Extract each item as a PyCell<PyFeature>
            let feature: &PyCell<PyFeature> = item.extract().expect("Failed to extract PyFeature");
            // Clone the inner Feature and keep the original index
            (index, feature.borrow().inner.clone())
        })
        .collect();

    // Set discriminant score to hyper score
    indexed_inner_collection.par_iter_mut().for_each(|(_, feat)| {

        match use_hyper_score {
            false => {
                feat.discriminant_score = (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0
            }
            true => {
                feat.discriminant_score = feat.hyperscore as f32;
            }
        }
    });

    // Sort indexed_inner_collection by discriminant_score
    indexed_inner_collection.par_sort_unstable_by(|(_, a), (_, b)| b.discriminant_score.total_cmp(&a.discriminant_score));

    // Extract the sorted indices
    let sorted_indices: Vec<usize> = indexed_inner_collection.iter().map(|(index, _)| *index).collect();

    // Perform additional operations on the sorted inner_collection
    let mut inner_collection: Vec<Feature> = indexed_inner_collection.into_iter().map(|(_, feat)| feat).collect();
    let _ = sage_core::ml::qvalue::spectrum_q_value(&mut inner_collection);
    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    // Update the original feature_collection according to the sorted order
    for (sorted_index, sorted_feature) in sorted_indices.iter().zip(inner_collection.iter()) {
        let feature: &PyCell<PyFeature> = feature_collection.get_item(*sorted_index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();
        // Update the feature's fields
        feature_borrow.inner.discriminant_score = sorted_feature.discriminant_score;
        feature_borrow.inner.spectrum_q = sorted_feature.spectrum_q;
        feature_borrow.inner.peptide_q = sorted_feature.peptide_q;
        feature_borrow.inner.protein_q = sorted_feature.protein_q;
    }

    Ok(())
}

#[pyfunction]
pub fn py_sage_fdr_psm(_py: Python, psm_collection: &PyList, indexed_database: &PyIndexedDatabase, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut indexed_inner_collection: Vec<(usize, Psm)> = psm_collection.iter()
        .enumerate()
        .map(|(index, item)| {
            // Extract each item as a PyCell<PyFeature>
            let feature: &PyCell<PyPsm> = item.extract().expect("Failed to extract PyFeature");
            // Clone the inner Feature and keep the original index
            (index, feature.borrow().inner.clone())
        })
        .collect();

    // Set discriminant score to hyper score
    indexed_inner_collection.par_iter_mut().for_each(|(_, feat)| {

        match use_hyper_score {
            false => {
                feat.sage_feature.discriminant_score = (-feat.sage_feature.poisson as f32).ln_1p() + feat.sage_feature.longest_y_pct / 3.0
            }
            true => {
                feat.sage_feature.discriminant_score = feat.sage_feature.hyperscore as f32;
            }
        }
    });

    // Sort indexed_inner_collection by discriminant_score
    indexed_inner_collection.par_sort_unstable_by(|(_, a), (_, b)| b.sage_feature.discriminant_score.total_cmp(&a.sage_feature.discriminant_score));

    // Extract the sorted indices
    let sorted_indices: Vec<usize> = indexed_inner_collection.iter().map(|(index, _)| *index).collect();

    // Perform additional operations on the sorted inner_collection
    let mut inner_collection: Vec<Feature> = indexed_inner_collection.into_iter().map(|(_, feat)| feat.sage_feature).collect();
    let _ = sage_core::ml::qvalue::spectrum_q_value(&mut inner_collection);
    let _ = picked_peptide(&indexed_database.inner, &mut inner_collection);
    let _ = picked_protein(&indexed_database.inner, &mut inner_collection);

    // Update the original feature_collection according to the sorted order
    for (sorted_index, sorted_feature) in sorted_indices.iter().zip(inner_collection.iter()) {
        let feature: &PyCell<PyPsm> = psm_collection.get_item(*sorted_index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();
        // Update the feature's fields
        feature_borrow.inner.sage_feature.discriminant_score = sorted_feature.discriminant_score;
        feature_borrow.inner.sage_feature.spectrum_q = sorted_feature.spectrum_q;
        feature_borrow.inner.sage_feature.peptide_q = sorted_feature.peptide_q;
        feature_borrow.inner.sage_feature.protein_q = sorted_feature.protein_q;
    }

    Ok(())
}

#[pymodule]
pub fn fdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    m.add_function(wrap_pyfunction!(py_sage_fdr, m)?)?;
    m.add_function(wrap_pyfunction!(py_sage_fdr_psm, m)?)?;
    Ok(())
}