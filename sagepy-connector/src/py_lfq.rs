use pyo3::prelude::*;
use sage_core::lfq::{IntegrationStrategy, PeakScoringStrategy};

#[pyclass]
pub struct PyPeakScoringStrategy {
    pub inner: PeakScoringStrategy,
}

#[pyclass]
pub struct PyIntegrationStrategy {
    pub inner: IntegrationStrategy,
}

#[pyclass]
pub struct P

#[pymodule]
pub fn lfq(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_class::<PyPeptideIx>()?;
    Ok(())
}
