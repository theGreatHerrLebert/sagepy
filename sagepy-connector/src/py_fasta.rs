use sage_core::fasta::Fasta;

use crate::py_enzyme::{PyDigest, PyEnzymeParameters};
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub struct PyFasta {
    pub inner: Fasta,
}

#[pymethods]
impl PyFasta {
    #[staticmethod]
    fn parse(contents: String, decoy_tag: String, generate_decoys: bool) -> PyResult<Self> {
        Ok(PyFasta {
            inner: Fasta::parse(contents, decoy_tag, generate_decoys),
        })
    }

    fn digest(&self, py: Python, enzyme_params: &PyEnzymeParameters) -> PyResult<PyObject> {
        let digests = self.inner.digest(&enzyme_params.inner);
        let py_digests: Vec<PyDigest> =
            digests.into_iter().map(|d| PyDigest { inner: d }).collect();
        Ok(py_digests.into_py(py))
    }
}

#[pymodule]
pub fn fasta(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFasta>()?;
    Ok(())
}
