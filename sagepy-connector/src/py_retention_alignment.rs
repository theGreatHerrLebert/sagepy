use pyo3::prelude::*;
use pyo3::types::PyList;
use sage_core::ml::retention_alignment::{Alignment, global_alignment};

use sage_core::scoring::Feature;
use crate::py_scoring::{PyFeature, PyPsm};

#[pyclass]
#[derive(Clone)]
pub struct PyAlignment {
    pub inner: Alignment,
}

#[pymethods]
impl PyAlignment {
    #[new]
    pub fn new(
        file_id: usize,
        max_rt: f32,
        slope: f32,
        intercept: f32,
    ) -> Self {
        PyAlignment {
            inner: Alignment {
                file_id,
                max_rt,
                slope,
                intercept,
            },
        }
    }
    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }
    #[getter]
    pub fn max_rt(&self) -> f32 {
        self.inner.max_rt
    }
    #[getter]
    pub fn slope(&self) -> f32 {
        self.inner.slope
    }
    #[getter]
    pub fn intercept(&self) -> f32 {
        self.inner.intercept
    }
}

#[pyfunction]
pub fn py_global_alignment(
    features: &Bound<'_, PyList>,
    n_files: usize,
) -> Vec<PyAlignment> {

    let mut inner_features: Vec<Feature> = features.iter()
        .map(|item| {
            let feature: Bound<'_, PyFeature> = item.extract().expect("Failed to extract PyFeature");
            feature.borrow().inner.clone()
        })
        .collect();

    global_alignment(&mut inner_features, n_files)
        .into_iter()
        .map(|alignment| PyAlignment { inner: alignment })
        .collect()
}

#[pyfunction]
pub fn py_global_alignment_psm(
    psms: &Bound<'_, PyList>,
    n_files: usize,
) -> Vec<PyAlignment> {
    // Step 1: clone out features + remember original index
    let indexed_psms: Vec<(usize, Feature)> = psms.iter()
        .enumerate()
        .map(|(i, item)| {
            let psm: Bound<'_, PyPsm> = item.extract().expect("Failed to extract PyPsm");
            (i, psm.borrow().inner.sage_feature.clone())
        })
        .collect();

    // Step 2: collect features to pass into alignment
    let mut features: Vec<Feature> = indexed_psms.iter().map(|(_, feat)| feat.clone()).collect();

    // Step 3: run global alignment on copied data
    let alignments = global_alignment(&mut features, n_files);

    // Step 4: write aligned_rt back into PyPsm
    for ((i, _), updated_feat) in indexed_psms.iter().zip(features.iter()) {
        let psm: Bound<'_, PyPsm> = psms.get_item(*i)
            .expect("Index out of range")
            .extract()
            .expect("Failed to extract PyPsm");
        let mut psm_mut = psm.borrow_mut();
        psm_mut.inner.sage_feature.aligned_rt = updated_feat.aligned_rt;
    }

    // Step 5: return the alignment parameters
    alignments
        .into_iter()
        .map(|alignment| PyAlignment { inner: alignment })
        .collect()
}

#[pymodule]
pub fn py_retention_alignment(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlignment>()?;
    m.add_function(wrap_pyfunction!(py_global_alignment, m)?)?;
    m.add_function(wrap_pyfunction!(py_global_alignment_psm, m)?)?;
    Ok(())
}