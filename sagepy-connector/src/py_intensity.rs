use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::intensity::FragmentIntensityPrediction;
use sage_core::intensity_prediction::{
    compute_key_hash, IntensityStore, PredictedIntensityStore, PredictedIntensityStoreV2,
};
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

/// Helper to convert u8 ion kind code to Kind enum.
fn parse_ion_kind(code: u8) -> Option<Kind> {
    match code {
        0 => Some(Kind::A),
        1 => Some(Kind::B),
        2 => Some(Kind::C),
        3 => Some(Kind::X),
        4 => Some(Kind::Y),
        5 => Some(Kind::Z),
        _ => None,
    }
}

/// Helper to parse ion kinds vector.
fn parse_ion_kinds(ion_kinds: Option<Vec<u8>>) -> Vec<Kind> {
    ion_kinds
        .unwrap_or_else(|| vec![1, 4]) // B and Y by default
        .iter()
        .filter_map(|&k| parse_ion_kind(k))
        .collect()
}

/// Python wrapper for the unified IntensityStore.
///
/// This class loads and provides access to pre-computed fragment ion intensity
/// predictions stored in a .sagi binary file. Supports both V1 (positional indexing)
/// and V2 (key-based indexing with sequence+charge keys) formats.
#[pyclass]
pub struct PyPredictedIntensityStore {
    pub inner: IntensityStore,
}

#[pymethods]
impl PyPredictedIntensityStore {
    /// Load predicted intensities from a .sagi binary file.
    ///
    /// Auto-detects the file format version (V1 or V2) and loads accordingly.
    #[staticmethod]
    pub fn load(path: &str) -> PyResult<Self> {
        let inner = IntensityStore::load(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Create a uniform V1 intensity store where all intensities are 1.0.
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
        let ion_kinds_parsed = parse_ion_kinds(ion_kinds);
        let inner = PredictedIntensityStore::uniform(&peptide_lengths, max_charge, ion_kinds_parsed);
        Self {
            inner: IntensityStore::V1(inner),
        }
    }

    /// Create a V2 intensity store from raw predictions.
    ///
    /// V2 format uses (sequence, charge) as the key, allowing database chunking
    /// and reuse of predictions across different database configurations.
    ///
    /// # Arguments
    /// * `entries` - List of tuples: (sequence_bytes, precursor_charge, peptide_len, intensities)
    ///   where intensities is a flat list of f32 values in [ion_kind][position][frag_charge] order
    /// * `max_charge` - Maximum fragment charge state
    /// * `ion_kinds` - Ion type codes (default: [1, 4] for B and Y)
    #[staticmethod]
    #[pyo3(signature = (entries, max_charge, ion_kinds=None))]
    pub fn from_raw_predictions_v2(
        entries: Vec<(Vec<u8>, u8, u16, Vec<f32>)>,
        max_charge: u8,
        ion_kinds: Option<Vec<u8>>,
    ) -> Self {
        let ion_kinds_parsed = parse_ion_kinds(ion_kinds);
        let inner = PredictedIntensityStoreV2::from_raw_predictions(
            entries,
            max_charge,
            ion_kinds_parsed,
        );
        Self {
            inner: IntensityStore::V2(inner),
        }
    }

    /// Compute the hash key for a (sequence, charge) pair.
    ///
    /// This uses the same hashing algorithm as the Rust implementation,
    /// ensuring consistency between Python and Rust lookups.
    #[staticmethod]
    pub fn compute_key_hash(sequence: &str, charge: u8) -> u64 {
        compute_key_hash(sequence.as_bytes(), charge)
    }

    /// Check if this is a V2 (key-based) store.
    #[getter]
    pub fn is_key_based(&self) -> bool {
        self.inner.is_key_based()
    }

    /// Get predicted intensity for a specific fragment (V1 positional lookup).
    ///
    /// # Arguments
    /// * `peptide_idx` - Index into IndexedDatabase.peptides
    /// * `peptide_len` - Length of the peptide sequence
    /// * `ion_kind` - Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
    /// * `position` - Fragment position (0-indexed)
    /// * `charge` - Charge state (1-indexed)
    ///
    /// # Returns
    /// The predicted intensity, or None if not found or if this is a V2 store.
    pub fn get_intensity(
        &self,
        peptide_idx: usize,
        peptide_len: usize,
        ion_kind: u8,
        position: usize,
        charge: u8,
    ) -> Option<f32> {
        let kind = parse_ion_kind(ion_kind)?;
        self.inner
            .get_intensity_by_idx(peptide_idx, peptide_len, kind, position, charge)
    }

    /// Get predicted intensity with default fallback (V1 positional lookup).
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
        self.get_intensity(peptide_idx, peptide_len, ion_kind, position, charge)
            .unwrap_or(1.0)
    }

    /// Get predicted intensity by sequence and charge (V2 key-based lookup).
    ///
    /// # Arguments
    /// * `sequence` - Raw peptide sequence (NOT UNIMOD notation)
    /// * `precursor_charge` - Precursor charge state
    /// * `ion_kind` - Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
    /// * `position` - Fragment position (0-indexed)
    /// * `fragment_charge` - Fragment charge state (1-indexed)
    ///
    /// # Returns
    /// The predicted intensity, or None if not found or if this is a V1 store.
    pub fn get_intensity_by_key(
        &self,
        sequence: &str,
        precursor_charge: u8,
        ion_kind: u8,
        position: usize,
        fragment_charge: u8,
    ) -> Option<f32> {
        let kind = parse_ion_kind(ion_kind)?;
        self.inner.get_intensity_by_key(
            sequence.as_bytes(),
            precursor_charge,
            kind,
            position,
            fragment_charge,
        )
    }

    /// Get predicted intensity by key with default fallback (V2 lookup).
    ///
    /// Returns 1.0 if the intensity is not found, which acts as a neutral weight.
    pub fn get_intensity_by_key_or_default(
        &self,
        sequence: &str,
        precursor_charge: u8,
        ion_kind: u8,
        position: usize,
        fragment_charge: u8,
    ) -> f32 {
        self.get_intensity_by_key(sequence, precursor_charge, ion_kind, position, fragment_charge)
            .unwrap_or(1.0)
    }

    /// Check if a key exists in the V2 store.
    ///
    /// Returns False for V1 stores.
    pub fn contains_key(&self, sequence: &str, precursor_charge: u8) -> bool {
        match &self.inner {
            IntensityStore::V2(store) => store.contains_key(sequence.as_bytes(), precursor_charge),
            IntensityStore::V1(_) => false,
        }
    }

    /// Get peptide length for a key (V2 only).
    ///
    /// Returns None for V1 stores or if key not found.
    pub fn get_peptide_len_by_key(&self, sequence: &str, precursor_charge: u8) -> Option<u16> {
        match &self.inner {
            IntensityStore::V2(store) => {
                store.get_peptide_len(sequence.as_bytes(), precursor_charge)
            }
            IntensityStore::V1(_) => None,
        }
    }

    /// Number of entries in the store.
    ///
    /// For V1, this is the number of peptides. For V2, this is the number of
    /// (sequence, charge) entries.
    #[getter]
    pub fn entry_count(&self) -> usize {
        match &self.inner {
            IntensityStore::V1(store) => store.peptide_count(),
            IntensityStore::V2(store) => store.entry_count(),
        }
    }

    /// Number of peptides in the store (V1 only, same as entry_count for V1).
    ///
    /// For V2 stores, this returns 0 (use entry_count instead).
    #[getter]
    pub fn peptide_count(&self) -> usize {
        match &self.inner {
            IntensityStore::V1(store) => store.peptide_count(),
            IntensityStore::V2(_) => 0, // V2 doesn't have peptide indices
        }
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

    /// Write the store to a .sagi binary file.
    ///
    /// The file format version (V1 or V2) is determined by the store type.
    pub fn write(&self, path: &str) -> PyResult<()> {
        match &self.inner {
            IntensityStore::V1(store) => store
                .write(path)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string())),
            IntensityStore::V2(store) => store
                .write(path)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string())),
        }
    }
}

#[pymodule]
pub fn py_intensity(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFragmentIntensityPrediction>()?;
    m.add_class::<PyPredictedIntensityStore>()?;
    Ok(())
}