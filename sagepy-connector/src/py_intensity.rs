use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::intensity::FragmentIntensityPrediction;
use crate::py_scoring::PyFragments;

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyFragmentIntensityPrediction {
    pub inner: FragmentIntensityPrediction,
}

#[pymethods]
impl PyFragmentIntensityPrediction {
    #[new]
    fn new(
        fragments: PyFragments,
        prosit_intensity_predicted: Vec<f32>,
    ) -> Self {
        PyFragmentIntensityPrediction {
            inner: FragmentIntensityPrediction {
                fragments: fragments.inner.clone(),
                prosit_intensity_predicted,
            },
        }
    }

    #[getter]
    fn prosit_intensity_predicted(&self) -> Vec<f32> {
        self.inner.prosit_intensity_predicted.clone()
    }

    #[setter]
    fn set_prosit_intensity_predicted(&mut self, prosit_intensity_predicted: Vec<f32>) {
        self.inner.prosit_intensity_predicted = prosit_intensity_predicted;
    }

    fn cosine_similarity(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        self.inner.cosine_similarity(epsilon, reduce_matched).unwrap()
    }

    fn spectral_angle_similarity(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        self.inner.spectral_angle_similarity(epsilon, reduce_matched)
    }

    fn pearson_correlation(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        self.inner.pearson_correlation(epsilon, reduce_matched)
    }

    fn spearman_correlation(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        self.inner.spearman_correlation(epsilon, reduce_matched)
    }

    fn spectral_entropy_similarity(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        self.inner.spectral_entropy_similarity(epsilon, reduce_matched)
    }

    fn observed_intensity_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        self.inner.observed_intensity_to_fragments_map()
    }

    fn predicted_intensity_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        self.inner.prosit_intensity_to_fragments_map()
    }

    fn prosit_intensity_to_fragments(&self) -> PyFragments {
        PyFragments {
            inner: self.inner.prosit_intensity_to_fragments(),
        }
    }
}

#[pymodule]
pub fn intensity(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragmentIntensityPrediction>()?;
    Ok(())
}