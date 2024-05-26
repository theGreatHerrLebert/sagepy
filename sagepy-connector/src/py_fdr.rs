use pyo3::prelude::*;
use sage_core::fdr::{Competition};
use sage_core::database::PeptideIx;
use crate::py_database::{PyPeptideIx};

#[pyclass]
// TODO: Check if it makes sense to tie this to PeptideIx
struct PyCompetitionPeptideIx {
    inner: Competition<PeptideIx>,
}

#[pymethods]
impl PyCompetitionPeptideIx {
    #[new]
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

#[pymodule]
pub fn fdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCompetitionPeptideIx>()?;
    Ok(())
}
