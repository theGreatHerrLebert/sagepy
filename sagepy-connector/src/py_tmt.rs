use pyo3::prelude::*;


#[pymodule]
pub fn tmt(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_class::<PyPeptideIx>()?;
    Ok(())
}
