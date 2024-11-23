use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::intensity::FragmentIntensityPrediction;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::py_scoring::PyPsm;

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyFragmentIntensityPrediction {
    pub inner: FragmentIntensityPrediction,
}

#[pymethods]
impl PyFragmentIntensityPrediction {
    #[new]
    fn new(
        intensities_observed: Vec<f32>,
        mz_observed: Vec<f32>,
        mz_calculated: Vec<f32>,
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
    fn mz_observed(&self) -> Vec<f32> {
        self.inner.mz_observed.clone()
    }

    #[setter]
    fn set_mz_observed(&mut self, mz_observed: Vec<f32>) {
        self.inner.mz_observed = mz_observed;
    }

    #[getter]
    fn mz_calculated(&self) -> Vec<f32> {
        self.inner.mz_calculated.clone()
    }

    #[setter]
    fn set_mz_calculated(&mut self, mz_calculated: Vec<f32>) {
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
}

#[pyfunction]
pub fn peptide_spectrum_match_to_feature_vector(
    psm: &PyPsm,
    epsilon: f32,
    reduce_matched: bool,
) -> Vec<f32> {
    let fragment_intensity_prediction = psm.inner.get_fragment_intensity_prediction();
    fragment_intensity_prediction.get_feature_vector(epsilon, reduce_matched)
}

#[pyfunction]
pub fn peptide_spectrum_match_list_to_intensity_feature_matrix_parallel(
    psms: Vec<PyPsm>,
    epsilon: f32,
    reduce_matched: bool,
    num_threads: usize,
) -> Vec<Vec<f32>> {
    let thread_pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    thread_pool.install(|| {
        psms.par_iter().map(|psm| {
            peptide_spectrum_match_to_feature_vector(psm, epsilon, reduce_matched)
        }).collect()
    })
}

#[pymodule]
pub fn intensity(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFragmentIntensityPrediction>()?;
    m.add_function(wrap_pyfunction!(peptide_spectrum_match_to_feature_vector, m)?)?;
    m.add_function(wrap_pyfunction!(peptide_spectrum_match_list_to_intensity_feature_matrix_parallel, m)?)?;
    Ok(())
}