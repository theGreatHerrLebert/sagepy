use pyo3::prelude::*;
use sage_core::tmt::Isobaric;

#[pyclass]
pub struct PyIsobaric {
    pub inner: Isobaric,
}

#[pymodule]
pub fn tmt(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIsobaric>()?;
    Ok(())
}
