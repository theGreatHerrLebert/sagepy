use std::collections::HashMap;
use pyo3::prelude::*;
use unimod::unimod::{unimod_modifications_mass_numerical, unimod_modifications_mass};

#[pyfunction]
fn unimod_modification_to_mass_numerical() -> HashMap<u32, f64> {
    unimod_modifications_mass_numerical()
}

#[pyfunction]
fn unimod_modification_to_mass() -> HashMap<&'static str, f64> {
    unimod_modifications_mass()
}

#[pymodule]
pub fn unimodifications(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass_numerical, m)?)?;
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass, m)?)?;
    Ok(())
}