use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::intensity::FragmentIntensityPrediction;
use sage_core::intensity_prediction::PredictedIntensityStore;
use sage_core::ion_series::Kind;
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

/// Python wrapper for the PredictedIntensityStore.
///
/// This class loads and provides access to pre-computed fragment ion intensity
/// predictions stored in a .sagi binary file.
#[pyclass]
pub struct PyPredictedIntensityStore {
    pub inner: PredictedIntensityStore,
}

#[pymethods]
impl PyPredictedIntensityStore {
    /// Load predicted intensities from a .sagi binary file.
    #[staticmethod]
    pub fn load(path: &str) -> PyResult<Self> {
        let inner = PredictedIntensityStore::load(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Create a uniform intensity store where all intensities are 1.0.
    ///
    /// This is useful for testing the weighted scoring code path without
    /// actual predictions. With uniform intensities, weighted scores should
    /// be identical to unweighted scores.
    ///
    /// # Arguments
    /// * `peptide_lengths` - Length of each peptide in the database
    /// * `max_charge` - Maximum fragment charge state (default: 3)
    /// * `ion_kinds` - Ion type codes (default: [1, 4] for B and Y)
    #[staticmethod]
    #[pyo3(signature = (peptide_lengths, max_charge=3, ion_kinds=None))]
    pub fn uniform(
        peptide_lengths: Vec<usize>,
        max_charge: u8,
        ion_kinds: Option<Vec<u8>>,
    ) -> Self {
        let ion_kinds_parsed: Vec<Kind> = ion_kinds
            .unwrap_or_else(|| vec![1, 4]) // B and Y by default
            .iter()
            .filter_map(|&k| match k {
                0 => Some(Kind::A),
                1 => Some(Kind::B),
                2 => Some(Kind::C),
                3 => Some(Kind::X),
                4 => Some(Kind::Y),
                5 => Some(Kind::Z),
                _ => None,
            })
            .collect();

        let inner = PredictedIntensityStore::uniform(&peptide_lengths, max_charge, ion_kinds_parsed);
        Self { inner }
    }

    /// Get predicted intensity for a specific fragment.
    ///
    /// # Arguments
    /// * `peptide_idx` - Index into IndexedDatabase.peptides
    /// * `peptide_len` - Length of the peptide sequence
    /// * `ion_kind` - Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
    /// * `position` - Fragment position (0-indexed)
    /// * `charge` - Charge state (1-indexed)
    ///
    /// # Returns
    /// The predicted intensity, or None if not found.
    pub fn get_intensity(
        &self,
        peptide_idx: usize,
        peptide_len: usize,
        ion_kind: u8,
        position: usize,
        charge: u8,
    ) -> Option<f32> {
        let kind = match ion_kind {
            0 => Kind::A,
            1 => Kind::B,
            2 => Kind::C,
            3 => Kind::X,
            4 => Kind::Y,
            5 => Kind::Z,
            _ => return None,
        };
        self.inner.get_intensity(peptide_idx, peptide_len, kind, position, charge)
    }

    /// Get predicted intensity with default fallback.
    ///
    /// Returns 1.0 if the intensity is not found, which acts as a neutral weight.
    pub fn get_intensity_or_default(
        &self,
        peptide_idx: usize,
        peptide_len: usize,
        ion_kind: u8,
        position: usize,
        charge: u8,
    ) -> f32 {
        let kind = match ion_kind {
            0 => Kind::A,
            1 => Kind::B,
            2 => Kind::C,
            3 => Kind::X,
            4 => Kind::Y,
            5 => Kind::Z,
            _ => return 1.0,
        };
        self.inner.get_intensity_or_default(peptide_idx, peptide_len, kind, position, charge)
    }

    /// Number of peptides in the store.
    #[getter]
    pub fn peptide_count(&self) -> usize {
        self.inner.peptide_count()
    }

    /// Maximum fragment charge state stored.
    #[getter]
    pub fn max_charge(&self) -> u8 {
        self.inner.max_charge()
    }

    /// Ion type codes stored (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z).
    #[getter]
    pub fn ion_kinds(&self) -> Vec<u8> {
        self.inner.ion_kinds().iter().map(|k| *k as u8).collect()
    }
}

#[pymodule]
pub fn py_intensity(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFragmentIntensityPrediction>()?;
    m.add_class::<PyPredictedIntensityStore>()?;
    Ok(())
}