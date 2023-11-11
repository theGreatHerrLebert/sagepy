use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use std::sync::Arc;

use std::hash::Hash;

use pyo3::exceptions::PyValueError;
use pyo3::types::PyList;
use sage_core::enzyme::{Digest, Enzyme, EnzymeParameters, Position};
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

#[pyclass]
#[derive(Clone)]
pub struct PyPosition {
    pub inner: Position,
}

#[pymethods]
impl PyPosition {
    #[staticmethod]
    fn nterm() -> Self {
        PyPosition {
            inner: Position::Nterm,
        }
    }

    #[staticmethod]
    fn cterm() -> Self {
        PyPosition {
            inner: Position::Cterm,
        }
    }

    #[staticmethod]
    fn full() -> Self {
        PyPosition {
            inner: Position::Full,
        }
    }

    #[staticmethod]
    fn internal() -> Self {
        PyPosition {
            inner: Position::Internal,
        }
    }

    #[staticmethod]
    fn from_string(position_string: &str) -> PyResult<Self> {
        match position_string {
            "n_term" => Ok(PyPosition::nterm()),
            "c_term" => Ok(PyPosition::cterm()),
            "full" => Ok(PyPosition::full()),
            "internal" => Ok(PyPosition::internal()),
            _ => Err(PyValueError::new_err("Invalid position string")),
        }
    }

    #[getter]
    fn to_string(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyDigest {
    pub inner: Digest,
}

#[pymethods]
impl PyDigest {
    #[new]
    fn new(
        decoy: bool,
        sequence: &str,
        protein: &str,
        missed_cleavages: u8,
        position: PyPosition,
        semi_enzymatic: bool,
    ) -> Self {
        PyDigest {
            inner: Digest {
                decoy,
                sequence: sequence.to_string(),
                protein: Arc::new(protein.to_string()),
                missed_cleavages,
                position: position.inner,
                semi_enzymatic,
            },
        }
    }

    #[getter]
    fn decoy(&self) -> bool {
        self.inner.decoy
    }

    #[getter]
    fn sequence(&self) -> &str {
        &self.inner.sequence
    }

    #[getter]
    fn protein(&self) -> &str {
        &self.inner.protein
    }

    #[getter]
    fn missed_cleavages(&self) -> u8 {
        self.inner.missed_cleavages
    }

    #[getter]
    fn position(&self) -> String {
        format!("{:?}", self.inner.position)
    }

    fn reverse(&self) -> PyResult<PyDigest> {
        Ok(PyDigest {
            inner: self.inner.reverse(),
        })
    }

    fn __eq__(&self, other: &PyDigest) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> isize {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish() as isize
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyEnzyme {
    pub inner: Enzyme,
}

#[pymethods]
impl PyEnzyme {
    #[new]
    fn new(
        cleave: &str,
        c_terminal: bool,
        semi_enzymatic: bool,
        skip_suffix: Option<char>,
    ) -> PyResult<Self> {
        match Enzyme::new(cleave, skip_suffix, c_terminal, semi_enzymatic) {
            Some(enzyme) => Ok(PyEnzyme { inner: enzyme }),
            None => Err(PyValueError::new_err("Failed to create Enzyme")),
        }
    }

    #[getter]
    fn c_terminal(&self) -> bool {
        self.inner.c_terminal
    }

    #[getter]
    fn skip_suffix(&self) -> Option<char> {
        self.inner.skip_suffix
    }

    #[getter]
    fn semi_enzymatic(&self) -> bool {
        self.inner.semi_enzymatic
    }

    fn cleavage_sites(&self, py: Python, sequence: &str) -> PyResult<Py<PyArray2<usize>>> {
        // Call the original cleavage_sites method
        let sites = self.inner.cleavage_sites(sequence);

        // Convert the Vec<Range<usize>> to Vec<usize> while flattening
        let sites_flat: Vec<usize> = sites
            .into_iter()
            .flat_map(|s| vec![s.site.start, s.site.end])
            .collect();

        let rows = sites_flat.len() / 2;
        let np_array: Py<PyArray2<usize>> =
            sites_flat.into_pyarray(py).reshape([rows, 2])?.to_owned();

        Ok(np_array)
    }
}

#[pyclass]
pub struct PyEnzymeParameters {
    pub inner: EnzymeParameters,
}

#[pymethods]
impl PyEnzymeParameters {
    #[new]
    fn new(missed_cleavages: u8, min_len: usize, max_len: usize, enzyme: Option<PyEnzyme>) -> Self {
        PyEnzymeParameters {
            inner: EnzymeParameters {
                missed_cleavages,
                min_len,
                max_len,
                enyzme: enzyme.map(|e| e.inner),
            },
        }
    }

    #[getter]
    fn missed_cleavages(&self) -> u8 {
        self.inner.missed_cleavages
    }

    #[getter]
    fn min_len(&self) -> usize {
        self.inner.min_len
    }

    #[getter]
    fn max_len(&self) -> usize {
        self.inner.max_len
    }

    #[getter]
    fn enzyme(&self, _py: Python) -> PyResult<Option<PyEnzyme>> {
        match &self.inner.enyzme {
            Some(enzyme) => Ok(Some(PyEnzyme {
                inner: enzyme.clone(),
            })),
            None => Ok(None),
        }
    }
    fn cleavage_sites(&self, py: Python, sequence: &str) -> PyResult<Py<PyArray2<usize>>> {
        // Call the original cleavage_sites method
        let sites = self.inner.cleavage_sites(sequence);

        // Convert the Vec<Range<usize>> to Vec<usize> while flattening
        let sites_flat: Vec<usize> = sites
            .into_iter()
            .flat_map(|s| vec![s.site.start, s.site.end])
            .collect();

        let rows = sites_flat.len() / 2;
        let np_array: Py<PyArray2<usize>> =
            sites_flat.into_pyarray(py).reshape([rows, 2])?.to_owned();

        Ok(np_array)
    }

    pub fn digest(&self, py: Python, sequence: &str, protein: &str) -> PyResult<Py<PyList>> {
        let digests = self.inner.digest(sequence, Arc::new(protein.to_string()));

        // Create an empty Python list
        let list: Py<PyList> = PyList::empty(py).into();

        // Iterate over the digests and append them to the list
        for digest in digests {
            let py_digest = Py::new(py, PyDigest { inner: digest })?;
            list.as_ref(py).append(py_digest)?;
        }

        Ok(list.into())
    }
}

#[pymodule]
pub fn enzyme(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyDigest>()?;
    m.add_class::<PyPosition>()?;
    m.add_class::<PyEnzyme>()?;
    m.add_class::<PyEnzymeParameters>()?;
    Ok(())
}
