use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

use crate::py_enzyme::{PyDigest, PyPosition};
use sage_core::peptide::Peptide;

#[pyclass]
#[derive(Clone)]
pub struct PyPeptide {
    pub inner: Peptide,
}

#[pymethods]
impl PyPeptide {
    #[new]
    pub fn new(
        decoy: bool,
        sequence: String,
        modifications: Vec<f32>,
        mono_isotopic: f32,
        missed_cleavages: u8,
        position: PyPosition,
        proteins: Vec<String>,
        semi_enzymatic: bool,
        n_term: Option<f32>,
        c_term: Option<f32>,
    ) -> PyResult<PyPeptide> {
        let sequence_bytes = sequence.into_bytes(); // Convert the string to Vec<u8>
        let boxed_sequence = sequence_bytes.into_boxed_slice(); // Convert Vec<u8> to Box<[u8]>
        let arc_sequence = Arc::from(boxed_sequence); // Convert Box<[u8]> to Arc<[u8]> without dereferencing

        // Convert Python list of strings to Vec<Arc<String>>
        let arc_proteins = proteins.into_iter().map(Arc::new).collect();

        Ok(PyPeptide {
            inner: Peptide {
                decoy,
                sequence: arc_sequence,
                modifications: modifications,
                nterm: n_term,
                cterm: c_term,
                monoisotopic: mono_isotopic,
                missed_cleavages,
                position: position.inner,
                proteins: arc_proteins,
                semi_enzymatic,
            },
        })
    }

    #[staticmethod]
    fn try_new_from_digest(digest: &PyDigest) -> PyResult<Self> {
        let peptide = Peptide::try_from(digest.inner.clone())
            .map_err(|_e| PyErr::new::<PyValueError, _>(format!("Error creating peptide.")))?;
        Ok(PyPeptide { inner: peptide })
    }

    #[getter]
    pub fn decoy(&self) -> bool {
        self.inner.decoy
    }

    #[getter]
    pub fn sequence(&self) -> &str {
        std::str::from_utf8(&self.inner.sequence).unwrap()
    }

    #[getter]
    pub fn modifications(&self) -> Vec<f32> {
        self.inner.modifications.clone()
    }

    #[getter]
    pub fn n_term(&self) -> Option<f32> {
        self.inner.nterm
    }

    #[getter]
    pub fn c_term(&self) -> Option<f32> {
        self.inner.cterm
    }

    #[getter]
    pub fn monoisotopic(&self) -> f32 {
        self.inner.monoisotopic
    }

    #[getter]
    pub fn missed_cleavages(&self) -> u8 {
        self.inner.missed_cleavages
    }

    #[getter]
    pub fn position(&self) -> PyPosition {
        PyPosition {
            inner: self.inner.position,
        }
    }

    #[getter]
    pub fn proteins(&self) -> Vec<String> {
        self.inner.proteins.iter().map(|s| s.to_string()).collect()
    }

    pub fn reverse(&self) -> PyPeptide {
        PyPeptide {
            inner: self.inner.reverse(),
        }
    }
}

#[pymodule]
pub fn peptide(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptide>()?;
    Ok(())
}
