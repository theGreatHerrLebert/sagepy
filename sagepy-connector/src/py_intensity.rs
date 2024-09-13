use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::intensity::FragmentIntensityPrediction;

#[pyclass]
pub struct PyFragmentIntensityPrediction {
    pub inner: FragmentIntensityPrediction,
}

#[pymethods]
impl PyFragmentIntensityPrediction {
    #[new]
    fn new(
        intensities_observed: Vec<f32>,
        mz_observed: Vec<f64>,
        mz_calculated: Vec<f64>,
        charges: Vec<i32>,
        ordinals: Vec<i32>,
        ion_types: Vec<bool>,
        prosit_intensity_predicted: Vec<f32>,
    ) -> Self {
        PyFragmentIntensityPrediction {
            inner: FragmentIntensityPrediction {
                intensities_observed,
                mz_observed,
                mz_calculated,
                charges,
                ordinals,
                ion_types,
                prosit_intensity_predicted,
            },
        }
    }

    #[getter]
    fn intensities_observed(&self) -> Vec<f32> {
        self.inner.intensities_observed.clone()
    }

    #[setter]
    fn set_intensities_observed(&mut self, intensities_observed: Vec<f32>) {
        self.inner.intensities_observed = intensities_observed;
    }

    #[getter]
    fn mz_observed(&self) -> Vec<f64> {
        self.inner.mz_observed.clone()
    }

    #[setter]
    fn set_mz_observed(&mut self, mz_observed: Vec<f64>) {
        self.inner.mz_observed = mz_observed;
    }

    #[getter]
    fn mz_calculated(&self) -> Vec<f64> {
        self.inner.mz_calculated.clone()
    }

    #[setter]
    fn set_mz_calculated(&mut self, mz_calculated: Vec<f64>) {
        self.inner.mz_calculated = mz_calculated;
    }

    #[getter]
    fn charges(&self) -> Vec<i32> {
        self.inner.charges.clone()
    }

    #[setter]
    fn set_charges(&mut self, charges: Vec<i32>) {
        self.inner.charges = charges;
    }

    #[getter]
    fn ordinals(&self) -> Vec<i32> {
        self.inner.ordinals.clone()
    }

    #[setter]
    fn set_ordinals(&mut self, ordinals: Vec<i32>) {
        self.inner.ordinals = ordinals;
    }

    #[getter]
    fn ion_types(&self) -> Vec<bool> {
        self.inner.ion_types.clone()
    }

    #[setter]
    fn set_ion_types(&mut self, ion_types: Vec<bool>) {
        self.inner.ion_types = ion_types;
    }

    #[getter]
    fn prosit_intensity_predicted(&self) -> Vec<f32> {
        self.inner.prosit_intensity_predicted.clone()
    }

    #[setter]
    fn set_prosit_intensity_predicted(&mut self, prosit_intensity_predicted: Vec<f32>) {
        self.inner.prosit_intensity_predicted = prosit_intensity_predicted;
    }

    fn cosine_similarity(&self, epsilon: f32) -> f32 {
        self.inner.cosine_similarity(epsilon).unwrap()
    }

    fn spectral_angle_similarity(&self, epsilon: f32) -> f32 {
        self.inner.spectral_angle_similarity(epsilon)
    }

    fn pearson_correlation(&self, epsilon: f32) -> f32 {
        self.inner.pearson_correlation(epsilon)
    }

    fn spearman_correlation(&self, epsilon: f32) -> f32 {
        self.inner.spearman_correlation(epsilon)
    }

    fn spectral_entropy_similarity(&self, epsilon: f32) -> f32 {
        self.inner.spectral_entropy_similarity(epsilon)
    }

    fn observed_intensity_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        self.inner.observed_intensity_to_fragments_map()
    }

    fn predicted_intensity_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        self.inner.prosit_intensity_to_fragments_map()
    }
}

#[pymodule]
pub fn intensity(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragmentIntensityPrediction>()?;
    Ok(())
}