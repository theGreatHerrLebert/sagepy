use pyo3::prelude::*;
use pyo3::types::PyList;
use qfdrust::psm::Psm;
use sage_core::fdr::{Competition, picked_peptide, picked_protein};

use sage_core::database::{PeptideIx};
use sage_core::scoring::Feature;
use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_mass::PyTolerance;
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
    #[pyo3(signature = (forward, reverse, forward_ix=None, reverse_ix=None))]
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
pub fn py_sage_fdr(_py: Python, feature_collection: &Bound<'_, PyList>, indexed_database: &PyIndexedDatabase, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut indexed_inner_collection: Vec<(usize, Feature)> = feature_collection.iter()
        .enumerate()
        .map(|(index, item)| {
            // Extract each item as a Bound<PyFeature>
            let feature: Bound<'_, PyFeature> = item.extract().expect("Failed to extract PyFeature");
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
        let feature: Bound<'_, PyFeature> = feature_collection.get_item(*sorted_index).expect("Failed to get PyFeature").extract()?;
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
pub fn py_sage_fdr_psm(_py: Python, psm_collection: &Bound<'_, PyList>, indexed_database: &PyIndexedDatabase, use_hyper_score: bool) -> PyResult<()> {

    // Extract the inner collection of Feature objects along with their original indices
    let mut indexed_inner_collection: Vec<(usize, Psm)> = psm_collection.iter()
        .enumerate()
        .map(|(index, item)| {
            // Extract each item as a Bound<PyPsm>
            let feature: Bound<'_, PyPsm> = item.extract().expect("Failed to extract PyFeature");
            // Clone the inner Feature and keep the original index
            (index, feature.borrow().inner.clone())
        })
        .collect();

    // Set discriminant score to hyper score
    indexed_inner_collection.par_iter_mut().for_each(|(_, feat)| {
        match use_hyper_score {
            false => {
                feat.sage_feature.discriminant_score = feat.re_score.unwrap_or(0.0) as f32;
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

    // Update the original psm_collection according to the sorted order
    for (sorted_index, sorted_feature) in sorted_indices.iter().zip(inner_collection.iter()) {
        let feature: Bound<'_, PyPsm> = psm_collection.get_item(*sorted_index).expect("Failed to get PyFeature").extract()?;
        let mut feature_borrow = feature.borrow_mut();
        feature_borrow.inner.sage_feature.discriminant_score = sorted_feature.discriminant_score;
        feature_borrow.inner.sage_feature.spectrum_q = sorted_feature.spectrum_q;
        feature_borrow.inner.sage_feature.peptide_q = sorted_feature.peptide_q;
        feature_borrow.inner.sage_feature.protein_q = sorted_feature.protein_q;
    }

    Ok(())
}

/// Fit sage's linear discriminant model on the provided PSMs and write the
/// resulting `discriminant_score` and `posterior_error` back onto each PSM.
/// Mirrors `sage_core::ml::linear_discriminant::score_psms`, which is the
/// step `sage` CLI runs when `predict_rt: true`. After this you can call
/// `py_sage_fdr_psm` with `use_hyper_score=false` to get discriminant-based
/// q-values matching the sage CLI.
#[pyfunction]
pub fn py_lda_score_psms(
    _py: Python,
    psm_collection: &Bound<'_, PyList>,
    precursor_tol: &PyTolerance,
) -> PyResult<bool> {
    let mut features: Vec<Feature> = psm_collection
        .iter()
        .map(|item| {
            let psm: Bound<'_, PyPsm> = item.extract().expect("Failed to extract PyPsm");
            psm.borrow().inner.sage_feature.clone()
        })
        .collect();

    let fitted = sage_core::ml::linear_discriminant::score_psms(
        &mut features,
        precursor_tol.inner,
    )
    .is_some();

    if !fitted {
        // Same heuristic fallback sage CLI uses when LDA fails.
        for feat in features.iter_mut() {
            feat.discriminant_score =
                (-feat.poisson as f32).ln_1p() + feat.longest_y_pct / 3.0;
        }
    }

    for (item, feat) in psm_collection.iter().zip(features.into_iter()) {
        let psm: Bound<'_, PyPsm> = item.extract()?;
        let mut borrow = psm.borrow_mut();
        borrow.inner.sage_feature.discriminant_score = feat.discriminant_score;
        borrow.inner.sage_feature.posterior_error = feat.posterior_error;
        // Also stash on `re_score` so py_sage_fdr_psm(use_hyper_score=false),
        // which routes through re_score, picks up the LDA discriminant instead
        // of overwriting it with 0.
        borrow.inner.re_score = Some(feat.discriminant_score as f64);
    }

    Ok(fitted)
}

/// Assign initial spectrum q-values by sorting PSMs on their raw Poisson score.
///
/// Sage uses this pre-fit q-value pass to choose the 1% target set that trains
/// the RT and mobility predictors when `predict_rt` is enabled.
#[pyfunction]
pub fn py_assign_initial_spectrum_q_psms(
    _py: Python,
    psm_collection: &Bound<'_, PyList>,
) -> PyResult<()> {
    let mut indexed_inner_collection: Vec<(usize, Psm)> = psm_collection
        .iter()
        .enumerate()
        .map(|(index, item)| {
            let psm: Bound<'_, PyPsm> = item.extract().expect("Failed to extract PyPsm");
            (index, psm.borrow().inner.clone())
        })
        .collect();

    indexed_inner_collection.par_sort_unstable_by(|(_, a), (_, b)| {
        a.sage_feature.poisson.total_cmp(&b.sage_feature.poisson)
    });

    let mut inner_collection: Vec<Feature> = indexed_inner_collection
        .iter()
        .map(|(_, feat)| feat.sage_feature.clone())
        .collect();
    let _ = sage_core::ml::qvalue::spectrum_q_value(&mut inner_collection);

    for (sorted_index, sorted_feature) in indexed_inner_collection
        .iter()
        .zip(inner_collection.iter())
    {
        let psm: Bound<'_, PyPsm> = psm_collection
            .get_item(sorted_index.0)
            .expect("Failed to get PyPsm")
            .extract()?;
        let mut feature_borrow = psm.borrow_mut();
        feature_borrow.inner.sage_feature.spectrum_q = sorted_feature.spectrum_q;
    }

    Ok(())
}

#[pymodule]
pub fn py_fdr(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    m.add_function(wrap_pyfunction!(py_sage_fdr, m)?)?;
    m.add_function(wrap_pyfunction!(py_sage_fdr_psm, m)?)?;
    m.add_function(wrap_pyfunction!(py_lda_score_psms, m)?)?;
    m.add_function(wrap_pyfunction!(py_assign_initial_spectrum_q_psms, m)?)?;
    Ok(())
}
